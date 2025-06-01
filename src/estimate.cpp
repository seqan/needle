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
                result.push_back(number);
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
        prev_counts.resize(size, std::vector<float>(num_experiments));
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

        if (!counters_initialised)
            counters_initialised = init_counter(seqs.size());

#pragma omp parallel for
        for (size_t i = 0; i < seqs.size(); ++i)
        {
            std::vector<uint64_t> const minimiser = get_minimiser(seqs[i], args);
            float const minimiser_count = minimiser.size();
            float const minimiser_pos = minimiser_count / 2.0f;

            auto agent = ibf.counting_agent();
            std::vector<float> all_counts(total_bins);
            std::ranges::copy(agent.bulk_count(minimiser), all_counts.begin());

            // For each experiment (file)
            for (size_t exp_idx = 0; exp_idx < num_experiments; ++exp_idx)
            {
                float prev_count = 0.0f;

                // Go through levels from highest threshold to lowest
                for (size_t l = 0; l < num_levels; ++l)
                {
                    size_t const level = num_levels - 1 - l; // Start from highest level
                    size_t const bin_idx = level * num_experiments + exp_idx;

                    // Check if this bin is deleted
                    if (std::ranges::find(deleted, bin_idx) != deleted.end())
                        continue;

                    // Get the count for this bin
                    float count = all_counts[bin_idx];

                    // False positive correction
                    float const expected_fp = minimiser_count * fprs[level][exp_idx];
                    float const corrected = count - expected_fp;
                    float const normalized =
                        std::max(0.0f, static_cast<float>(corrected / (1.0 - fprs[level][exp_idx])));

                    // Check if we've found enough minimizers to reach the median
                    if (prev_count + normalized >= minimiser_pos)
                    {
                        // Level found
                        if (l == num_levels - 1) // This is the last (lowest) level
                        {
                            if constexpr (samplewise)
                                estimations[i][exp_idx] = expressions[level][exp_idx];
                            else
                                estimations[i][exp_idx] = args.expression_thresholds[level];
                        }
                        else
                        {
                            // Interpolate between this level and the previous one
                            assert(minimiser_pos > prev_count);
                            assert(normalized > 0.0f);

                            float const normalized_minimiser_pos = (minimiser_pos - prev_count) / normalized;

                            if constexpr (samplewise)
                            {
                                size_t const prev_level = level + 1; // Previous level had higher threshold
                                uint16_t const current_threshold = expressions[level][exp_idx];
                                uint16_t const prev_threshold = expressions[prev_level][exp_idx];
                                uint16_t const threshold_diff = prev_threshold - current_threshold;
                                uint16_t const estimate =
                                    prev_threshold - static_cast<uint16_t>(normalized_minimiser_pos * threshold_diff);
                                estimations[i][exp_idx] = std::max(current_threshold, estimate);
                            }
                            else
                            {
                                uint16_t const current_threshold = args.expression_thresholds[level];
                                uint16_t const prev_threshold = args.expression_thresholds[level + 1];
                                uint16_t const threshold_diff = prev_threshold - current_threshold;
                                uint16_t const estimate =
                                    prev_threshold - static_cast<uint16_t>(normalized_minimiser_pos * threshold_diff);
                                estimations[i][exp_idx] = std::max(current_threshold, estimate);
                            }
                        }

                        // Apply normalization if requested
                        if constexpr (normalization_method && samplewise)
                            estimations[i][exp_idx] = static_cast<uint16_t>(
                                estimations[i][exp_idx] / expressions[0][exp_idx]); // Normalize by first level

                        break; // Found the estimate for this experiment
                    }
                    else
                    {
                        // Add to previous count and continue to next level
                        prev_count += normalized;
                    }
                }
            }
        }

        // Write results: one line per sequence, one value per experiment
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
