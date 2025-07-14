// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "estimate.hpp"

#include <omp.h>

#include <cereal/archives/binary.hpp>

#include "misc/debug.hpp"
#include "misc/filenames.hpp"
#include "misc/read_levels.hpp"

inline std::vector<uint64_t> get_minimiser(seqan3::dna4_vector const & seq, minimiser_arguments const & args)
{
    std::vector<uint64_t> minimiser;
    auto view = seqan3::views::minimiser_hash(seq, args.shape, args.w_size, args.s) | std::views::common;
    // Or `return std::ranges::to<std::vector<uint64_t>>(view);`
    minimiser.assign(view.begin(), view.end());
    return minimiser;
}

// Actual estimation
template <typename ibf_t, bool normalization, typename exp_t>
void check_ibf(estimate_ibf_arguments const & args,
               ibf_t const & ibf,
               std::vector<uint16_t> & estimations,
               seqan3::dna4_vector const & seq,
               exp_t const & expressions,
               std::vector<std::vector<double>> const & fprs,
               std::vector<uint64_t> const & deleted,
               size_t const num_levels,
               size_t const num_experiments)
{
    // Check, if one expression threshold for all or individual thresholds
    static constexpr bool multiple_expressions = std::same_as<exp_t, std::vector<std::vector<uint16_t>>>;

    std::vector<uint64_t> const minimiser = get_minimiser(seq, args);
    uint64_t const minimiser_count = minimiser.size();

#ifndef NDEBUG
    if (minimiser_count > std::numeric_limits<uint16_t>::max())
        log("Warning: 16-bit counter might overflow.");
#endif

    std::vector<float> counter(ibf.bin_count());
    auto agent = ibf.counting_agent();
    std::ranges::copy(agent.bulk_count(minimiser), counter.begin());

    // Defines where the median should be
    float const minimiser_pos = minimiser_count / 2.0;

    // For each experiment (file)
    for (size_t experiment = 0; experiment < num_experiments; ++experiment)
    {
        float prev_count = 0.0f;

        // Go through levels from highest threshold to lowest
        for (size_t const level : std::views::iota(0u, num_levels) | std::views::reverse)
        {
            size_t const bin = level * num_experiments + experiment;

            // Check if this bin is deleted
            if (std::ranges::find(deleted, bin) != deleted.end())
                continue;

            // Correction by substracting the expected number of false positives
            counter[bin] = [&]()
            {
                float const expected_false_positives = minimiser_count * fprs[level][experiment];
                float const corrected_count = counter[bin] - expected_false_positives;
                float const normalized_count = std::max(0.0, corrected_count / (1.0 - fprs[level][experiment]));
                return normalized_count;
            }();

            // Check if considering previously seen minimisers and minimisers found at current level equal to or are greater
            // than the minimiser_pow, which gives the median position.
            // If an estimation took already place (estimations[experiment]!=0), a second estimation is not performed.
            if (estimations[experiment] == 0 && prev_count + counter[bin] >= minimiser_pos)
            {
                // Level found
                if (level == num_levels - 1) // This is the last (lowest) level
                {
                    if constexpr (multiple_expressions)
                        estimations[experiment] = expressions[level][experiment];
                    else
                        estimations[experiment] = args.expression_thresholds[level];
                }
                else
                {
                    // Interpolate between this level and the previous one
                    assert(minimiser_pos > prev_count);
                    // Should never fail. This would mean that prev_counts[bin] was enough by itself and we should already
                    // have estimated on the previous level.
                    assert(counter[bin] > 0.0f);
                    float const normalized_minimiser_pos = (minimiser_pos - prev_count) / counter[bin];

                    // Actually calculate estimation, in the else case level stands for the prev_expression
                    if constexpr (multiple_expressions)
                    {
                        size_t const prev_level_expression = expressions[level + 1][experiment];
                        size_t const expression_difference = prev_level_expression - expressions[level][experiment];
                        size_t const estimate =
                            prev_level_expression - (normalized_minimiser_pos * expression_difference);
                        estimations[experiment] = std::max<size_t>(expressions[level][experiment], estimate);
                    }
                    else
                    {
                        size_t const prev_level_expression = args.expression_thresholds[level + 1];
                        size_t const expression_difference = prev_level_expression - args.expression_thresholds[level];
                        size_t const estimate =
                            prev_level_expression - (normalized_minimiser_pos * expression_difference);
                        estimations[experiment] = std::max<size_t>(args.expression_thresholds[level], estimate);
                    }
                }

                // Apply normalization if requested
                // TODO: Is this meant to be expressions[0]?
                if constexpr (normalization && multiple_expressions)
                    estimations[experiment] /= expressions[1][experiment]; // Normalize by first level

                break; // Found the estimate for this experiment
            }
            else
            {
                // Add to previous count and continue to next level
                prev_count += counter[bin];
            }
        }
    }
}

/*! \brief Function to estimate expression value.
*  \param args        The arguments.
*  \param ibf         The ibf determing what kind ibf is used (compressed or uncompressed).
*  \param file_out    The output file.
*  \param estimate_args  The estimate arguments.
*/
template <typename ibf_t, bool samplewise, bool normalization_method = false>
void estimate(estimate_ibf_arguments & args,
              ibf_t & ibf,
              std::filesystem::path const & file_out,
              estimate_arguments const & estimate_args)
{
    // ========================================================================
    // const data
    // ========================================================================
    std::vector<std::vector<uint16_t>> const expressions = [&]()
    {
        std::vector<std::vector<uint16_t>> result;
        if constexpr (samplewise)
            read_levels<uint16_t>(result, filenames::levels(estimate_args.path_in));
        return result;
    }();

    std::vector<std::vector<double>> const fprs = [&]()
    {
        std::vector<std::vector<double>> result;
        read_levels<double>(result, filenames::fprs(estimate_args.path_in));
        return result;
    }();

    std::vector<uint64_t> const deleted = [&]()
    {
        std::vector<uint64_t> result;
        if (std::filesystem::path const deleted_files_path = filenames::deleted(estimate_args.path_in);
            std::filesystem::exists(deleted_files_path))
        {
            std::ifstream fin{deleted_files_path};
            uint64_t number;

            while (fin >> number)
            {
                result.push_back(number);
            }
        }
        return result;
    }();

    // ========================================================================
    // data used for the estimation
    // ========================================================================
    std::vector<std::string> ids;
    std::vector<seqan3::dna4_vector> seqs;
    std::vector<std::vector<float>> prev_counts;
    std::vector<std::vector<uint16_t>> estimations;
    bool counters_initialised = false;

    // ========================================================================
    // I/O
    // ========================================================================
    std::ofstream outfile{file_out};
    seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::id, seqan3::field::seq>> fin{
        estimate_args.search_file};

    // ========================================================================
    // Misc
    // ========================================================================
    omp_set_num_threads(args.threads);
    seqan3::contrib::bgzf_thread_count = args.threads; // I/O and OpenMP are separate
    std::ranges::sort(args.expression_thresholds);

    // ========================================================================
    // Load the single IBF
    // ========================================================================
    load_ibf(ibf, filenames::ibf(estimate_args.path_in, samplewise, 0, args));

    // Calculate dimensions based on IBF structure
    size_t const num_levels = args.number_expression_thresholds;
    size_t const total_bins = ibf.bin_count();

    assert(total_bins % num_levels == 0);                   // Ensure that total_bins is divisible by num_levels
    size_t const num_experiments = total_bins / num_levels; // This should match the number of files/experiments

    // ========================================================================
    // Helper functions
    // ========================================================================
    auto init_counter = [&](size_t const size)
    {
        static_assert(std::same_as<std::ranges::range_value_t<decltype(prev_counts)>, std::vector<float>>);
        prev_counts.resize(size, std::vector<float>(num_experiments));

        static_assert(std::same_as<std::ranges::range_value_t<decltype(estimations)>, std::vector<uint16_t>>);
        estimations.resize(size, std::vector<uint16_t>(num_experiments));

        return true;
    };

    auto clear_data = [&]()
    {
        ids.clear();
        seqs.clear();
        std::ranges::for_each(prev_counts,
                              [](auto & v)
                              {
                                  std::ranges::fill(v, float{});
                              });
        std::ranges::for_each(estimations,
                              [](auto & v)
                              {
                                  std::ranges::fill(v, uint16_t{});
                              });
    };

    auto process_ibf = [&](size_t const i)
    {
        if constexpr (samplewise)
        {
            check_ibf<ibf_t, normalization_method>(args,
                                                   ibf,
                                                   estimations[i],
                                                   seqs[i],
                                                   expressions,
                                                   fprs,
                                                   deleted,
                                                   num_levels,
                                                   num_experiments);
        }
        else
        {
            check_ibf<ibf_t, false>(args,
                                    ibf,
                                    estimations[i],
                                    seqs[i],
                                    args.expression_thresholds,
                                    fprs,
                                    deleted,
                                    num_levels,
                                    num_experiments);
        }
    };

    // ========================================================================
    // Main algorithm
    // ========================================================================
    for (auto && chunked_records : fin | seqan::stl::views::chunk(estimate_args.batch_size))
    {
        clear_data();

        // Process chunk records
        for (auto & [id, seq] : chunked_records)
        {
            ids.push_back(id);
            seqs.push_back(seq);
        }

        // If there are less records than the batch size, only use as much as needed.
        if (!counters_initialised)
            counters_initialised = init_counter(seqs.size());

#pragma omp parallel for
        for (size_t i = 0; i < seqs.size(); ++i)
        {
            process_ibf(i);
        }

        // Write results
        for (size_t i = 0; i < seqs.size(); ++i)
        {
            outfile << ids[i] << '\t';
            for (size_t j = 0; j < num_experiments; ++j)
                outfile << estimations[i][j] << '\t';
            outfile << '\n';
        }
    }
}

// Calls the correct form of estimate
void call_estimate(estimate_ibf_arguments & args, estimate_arguments & estimate_args)
{
    load_args(args, filenames::data(estimate_args.path_in));

    auto call = [&]<typename ibf_t>(ibf_t && ibf)
    {
        if (args.samplewise)
        {
            if (estimate_args.normalization_method)
                estimate<ibf_t, true, true>(args, ibf, args.path_out, estimate_args);
            else
                estimate<ibf_t, true, false>(args, ibf, args.path_out, estimate_args);
        }
        else
        {
            estimate<ibf_t, false, false>(args, ibf, args.path_out, estimate_args);
        }
    };

    call.template operator()<seqan::hibf::interleaved_bloom_filter>({});
}
