// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "estimate.hpp"

#include <cereal/archives/binary.hpp>

#include "misc/filenames.hpp"
#include "misc/read_levels.hpp"

// Actual estimation
template <class IBFType, bool last_exp, bool normalization, typename exp_t>
void check_ibf(minimiser_arguments const & args,
               IBFType const & ibf,
               std::vector<uint16_t> & estimations_i,
               seqan3::dna4_vector const & seq,
               std::vector<float> & prev_counts,
               exp_t const & expressions,
               uint16_t const level,
               std::vector<double> const & fprs,
               std::vector<uint64_t> const & deleted)
{
    // Check, if one expression threshold for all or individual thresholds
    static constexpr bool multiple_expressions = std::same_as<exp_t, std::vector<std::vector<uint16_t>>>;

    // Count minimisers in ibf of current level
    std::vector<float> counter(ibf.bin_count());
    uint64_t minimiser_count{};
    auto agent = ibf.membership_agent();
    for (auto minHash : seqan3::views::minimiser_hash(seq, args.shape, args.w_size, args.s))
    {
        std::ranges::transform(counter, agent.bulk_contains(minHash), counter.begin(), std::plus<uint32_t>());
        ++minimiser_count;
    }

    // Defines where the median should be
    float const minimiser_pos = minimiser_count / 2.0;

    // Check every experiment by going over the number of bins in the ibf.
    for (size_t bin = 0; bin < counter.size(); ++bin)
    {
        if (std::ranges::find(deleted, bin) != deleted.end())
            continue;

        // Correction by substracting the expected number of false positives
        counter[bin] = [&]()
        {
            float const expected_false_positives = minimiser_count * fprs[bin];
            float const corrected_count = counter[bin] - expected_false_positives;
            float const normalized_count = std::max(0.0, corrected_count / (1.0 - fprs[bin]));
            return normalized_count;
        }();

        // Check if considering previously seen minimisers and minimisers found at current level equal to or are greater
        // than the minimiser_pow, which gives the median position.
        // If an estimation took already place (estimations_i[bin]!=0), a second estimation is not performed.
        if (estimations_i[bin] == 0 && prev_counts[bin] + counter[bin] >= minimiser_pos)
        {
            // If there was no previous level, because we are looking at the last level.
            if constexpr (last_exp)
            {
                if constexpr (multiple_expressions)
                    estimations_i[bin] = expressions[level][bin];
                else
                    estimations_i[bin] = expressions;
            }
            else
            {
                assert(minimiser_pos > prev_counts[bin]);
                // Should never fail. This would mean that prev_counts[bin] was enough by itself and we should already
                // have estimated on the previous level.
                assert(counter[bin] > 0.0);
                float const normalized_minimiser_pos = (minimiser_pos - prev_counts[bin]) / counter[bin];

                // Actually calculate estimation, in the else case level stands for the prev_expression
                if constexpr (multiple_expressions)
                {
                    size_t const prev_level_expression = expressions[level + 1][bin];
                    size_t const expression_difference = prev_level_expression - expressions[level][bin];
                    size_t const estimate = prev_level_expression - (normalized_minimiser_pos * expression_difference);

                    estimations_i[bin] = std::max<size_t>(expressions[level][bin], estimate);
                }
                else
                {
                    size_t const prev_level_expression = level;
                    size_t const expression_difference = prev_level_expression - expressions;
                    size_t const estimate = prev_level_expression - (normalized_minimiser_pos * expression_difference);

                    estimations_i[bin] = std::max<size_t>(expressions, estimate);
                }
            }

            // Perform normalization by dividing through the threshold of the first level. Only works if multiple expressions were used.
            if constexpr (normalization && multiple_expressions)
                estimations_i[bin] /= expressions[1][bin];
        }
        else
        {
            // If not found at this level, add to previous count.
            prev_counts[bin] += counter[bin];
        }
    }
}

/*! \brief Function to estimate expression value.
*  \param args        The arguments.
*  \param ibf         The ibf determing what kind ibf is used (compressed or uncompressed).
*  \param file_out    The output file.
*  \param estimate_args  The estimate arguments.
*/
template <class IBFType, bool samplewise, bool normalization_method = false>
void estimate(estimate_ibf_arguments & args,
              IBFType & ibf,
              std::filesystem::path file_out,
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
        if (std::filesystem::exists(filenames::deleted(estimate_args.path_in)))
        {
            std::ifstream fin{filenames::deleted(estimate_args.path_in)};
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
    uint64_t prev_expression{};
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
    // Helper functions
    // ========================================================================
    auto init_counter = [&](size_t const size)
    {
        static_assert(std::same_as<std::ranges::range_value_t<decltype(prev_counts)>, std::vector<float>>);
        prev_counts.resize(size, std::vector<float>(ibf.bin_count()));

        static_assert(std::same_as<std::ranges::range_value_t<decltype(estimations)>, std::vector<uint16_t>>);
        estimations.resize(size, std::vector<uint16_t>(ibf.bin_count()));

        return true;
    };

    auto clear_data = [&]()
    {
        prev_expression = 0u;
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

    auto process_ibf = [&]<bool is_last_level>(size_t const i, size_t const level)
    {
        if constexpr (samplewise)
        {
            check_ibf<IBFType, is_last_level, normalization_method>(args,
                                                                    ibf,
                                                                    estimations[i],
                                                                    seqs[i],
                                                                    prev_counts[i],
                                                                    expressions,
                                                                    level,
                                                                    fprs[level],
                                                                    deleted);
        }
        else
        {
            check_ibf<IBFType, is_last_level, false>(args,
                                                     ibf,
                                                     estimations[i],
                                                     seqs[i],
                                                     prev_counts[i],
                                                     args.expression_thresholds[level],
                                                     prev_expression,
                                                     fprs[level],
                                                     deleted);
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

        bool is_last_level = true;
        for (size_t const level : std::views::iota(0u, args.number_expression_thresholds) | std::views::reverse)
        {
            load_ibf(ibf, filenames::ibf(estimate_args.path_in, samplewise, level, args));

            // If there are less records than the batch size, only use as much as needed.
            if (!counters_initialised)
                counters_initialised = init_counter(seqs.size());

#pragma omp parallel for
            for (size_t i = 0; i < seqs.size(); ++i)
            {
                if (is_last_level)
                    process_ibf.template operator()<true>(i, level);
                else
                    process_ibf.template operator()<false>(i, level);
            }

            if constexpr (!samplewise)
                prev_expression = args.expression_thresholds[level];

            is_last_level = false;
        }

        // Write results
        for (size_t i = 0; i < seqs.size(); ++i)
        {
            outfile << ids[i] << '\t';
            for (size_t j = 0; j < ibf.bin_count(); ++j)
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

    if (args.compressed)
        call.template operator()<seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed>>({});
    else
        call.template operator()<seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>>({});
}
