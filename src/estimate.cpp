// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <algorithm>
#include <deque>
#include <iostream>
#include <math.h>
#include <numeric>
#include <omp.h>
#include <ranges>
#include <stdlib.h>
#include <string>
#include <vector>

#if SEQAN3_WITH_CEREAL
#    include <cereal/archives/binary.hpp>
#endif // SEQAN3_WITH_CEREAL

#include <seqan3/alphabet/container/concatenated_sequences.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/io/sequence_file/all.hpp>

#include "estimate.hpp"

// Actual estimation
template <class IBFType, bool last_exp, bool normalization, typename exp_t>
void check_ibf(min_arguments const & args,
               IBFType const & ibf,
               std::vector<uint16_t> & estimations_i,
               seqan3::dna4_vector const & seq,
               std::vector<uint32_t> & prev_counts,
               exp_t const & expressions,
               uint16_t const level,
               std::vector<double> const fprs,
               std::vector<int> & deleted)
{
    // Check, if one expression threshold for all or individual thresholds
    static constexpr bool multiple_expressions = std::same_as<exp_t, std::vector<std::vector<uint16_t>>>;

    // Count minimisers in ibf of current level
    std::vector<uint32_t> counter(ibf.bin_count());
    uint64_t minimiser_count{};
    auto agent = ibf.membership_agent();
    for (auto minHash : seqan3::views::minimiser_hash(seq, args.shape, args.w_size, args.s))
    {
        std::ranges::transform(counter, agent.bulk_contains(minHash), counter.begin(), std::plus<uint32_t>());
        ++minimiser_count;
    }

    // Defines where the median should be
    double const minimiser_pos = minimiser_count / 2.0;

    // Check every experiment by going over the number of bins in the ibf.
    for (size_t bin = 0; bin < counter.size(); ++bin)
    {
        if (std::ranges::find(deleted, bin) != deleted.end())
            continue;

        // Correction by substracting the expected number of false positives
        double const expected_false_positives = minimiser_count * fprs[bin];
        double const corrected_count = counter[bin] - expected_false_positives;
        double const normalized_count = corrected_count / (1.0 - fprs[bin]);

        // TODO: Rounding?
        // `normalized_count` could be ceiled
        // `estimate` further down could be ceiled
        counter[bin] = normalized_count < 0.0 ? 0u : normalized_count;

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
                assert(counter[bin] > 0u);
                double const normalized_minimiser_pos = (minimiser_pos - prev_counts[bin]) / counter[bin];

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

// Reads the level file ibf creates
// Defined in ibf.cpp
template <typename float_or_int>
void read_levels(std::vector<std::vector<float_or_int>> &, std::filesystem::path);

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
    std::vector<std::string> ids;
    std::vector<seqan3::dna4_vector> seqs;
    std::vector<uint32_t> counter;
    std::vector<uint16_t> counter_est;
    std::vector<std::vector<uint32_t>> prev_counts;
    uint64_t prev_expression;
    std::vector<std::vector<uint16_t>> estimations;
    std::vector<std::vector<uint16_t>> expressions;
    std::vector<std::vector<double>> fprs;
    std::vector<int> deleted{};

    omp_set_num_threads(args.threads);
    seqan3::contrib::bgzf_thread_count = args.threads;

    seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::id, seqan3::field::seq>> fin{
        estimate_args.search_file};
    for (auto & [id, seq] : fin)
    {
        ids.push_back(id);
        seqs.push_back(seq);
    }

    if constexpr (samplewise)
        read_levels<uint16_t>(expressions, estimate_args.path_in.string() + "IBF_Levels.levels");
    else
        prev_expression = 0;

    read_levels<double>(fprs, estimate_args.path_in.string() + "IBF_FPRs.fprs");

    // Check, if bins have been deleted
    if (std::filesystem::exists(estimate_args.path_in.string() + "IBF_Deleted"))
    {
        std::ifstream fin;
        fin.open(estimate_args.path_in.string() + "IBF_Deleted");
        uint64_t number;

        while (fin >> number)
        {
            deleted.push_back(number);
        }
        fin.close();
    }

    // Make sure expression levels are sorted.
    sort(args.expression_thresholds.begin(), args.expression_thresholds.end());

    // Initialse last expression.
    if constexpr (samplewise)
        load_ibf(ibf,
                 estimate_args.path_in.string() + "IBF_Level_" + std::to_string(args.number_expression_thresholds - 1));
    else
        load_ibf(ibf,
                 estimate_args.path_in.string() + "IBF_"
                     + std::to_string(args.expression_thresholds[args.expression_thresholds.size() - 1]));
    counter.assign(ibf.bin_count(), 0);
    counter_est.assign(ibf.bin_count(), 0);

    for (size_t i = 0; i < seqs.size(); ++i)
    {
        estimations.push_back(counter_est);
        prev_counts.push_back(counter);
    }
    counter_est.clear();
    counter.clear();

// Go over the sequences
#pragma omp parallel for
    for (size_t i = 0; i < seqs.size(); ++i)
    {
        if constexpr (samplewise && normalization_method)
            check_ibf<IBFType, true, true>(args,
                                           ibf,
                                           estimations[i],
                                           seqs[i],
                                           prev_counts[i],
                                           expressions,
                                           args.number_expression_thresholds - 1,
                                           fprs[args.number_expression_thresholds - 1],
                                           deleted);
        else if constexpr (samplewise)
            check_ibf<IBFType, true, false>(args,
                                            ibf,
                                            estimations[i],
                                            seqs[i],
                                            prev_counts[i],
                                            expressions,
                                            args.number_expression_thresholds - 1,
                                            fprs[args.number_expression_thresholds - 1],
                                            deleted);
        else
            check_ibf<IBFType, true, false>(args,
                                            ibf,
                                            estimations[i],
                                            seqs[i],
                                            prev_counts[i],
                                            args.expression_thresholds[args.expression_thresholds.size() - 1],
                                            prev_expression,
                                            fprs[args.expression_thresholds.size() - 1],
                                            deleted);
    }

    if constexpr (!samplewise)
        prev_expression = args.expression_thresholds[args.expression_thresholds.size() - 1];

    for (int j = args.number_expression_thresholds - 2; j >= 0; j--)
    {
        // Loadthe next ibf that should be considered.
        if constexpr (samplewise)
            load_ibf(ibf, estimate_args.path_in.string() + "IBF_Level_" + std::to_string(j));
        else
            load_ibf(ibf, estimate_args.path_in.string() + "IBF_" + std::to_string(args.expression_thresholds[j]));

// Go over the sequences
#pragma omp parallel for
        for (size_t i = 0; i < seqs.size(); ++i)
        {
            if constexpr (samplewise && normalization_method)
                check_ibf<IBFType, false, true>(args,
                                                ibf,
                                                estimations[i],
                                                seqs[i],
                                                prev_counts[i],
                                                expressions,
                                                j,
                                                fprs[j],
                                                deleted);
            else if constexpr (samplewise)
                check_ibf<IBFType, false, false>(args,
                                                 ibf,
                                                 estimations[i],
                                                 seqs[i],
                                                 prev_counts[i],
                                                 expressions,
                                                 j,
                                                 fprs[j],
                                                 deleted);
            else
                check_ibf<IBFType, false, false>(args,
                                                 ibf,
                                                 estimations[i],
                                                 seqs[i],
                                                 prev_counts[i],
                                                 args.expression_thresholds[j],
                                                 prev_expression,
                                                 fprs[j],
                                                 deleted);
        }

        if (!samplewise)
            prev_expression = args.expression_thresholds[j];
    }

    // Write output file.
    std::ofstream outfile;
    outfile.open(std::string{file_out});
    for (size_t i = 0; i < seqs.size(); ++i)
    {
        outfile << ids[i] << "\t";
        for (size_t j = 0; j < ibf.bin_count(); ++j)
            outfile << estimations[i][j] << "\t";

        outfile << "\n";
    }
    outfile.close();
}

// Calls the correct form of estimate
void call_estimate(estimate_ibf_arguments & args, estimate_arguments & estimate_args)
{
    load_args(args, std::string{estimate_args.path_in} + "IBF_Data");

    if (args.compressed)
    {
        seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
        if (args.samplewise)
        {
            if (estimate_args.normalization_method)
                estimate<seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed>, true, true>(args,
                                                                                                        ibf,
                                                                                                        args.path_out,
                                                                                                        estimate_args);
            else
                estimate<seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed>, true>(args,
                                                                                                  ibf,
                                                                                                  args.path_out,
                                                                                                  estimate_args);
        }
        else
        {
            estimate<seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed>, false>(args,
                                                                                               ibf,
                                                                                               args.path_out,
                                                                                               estimate_args);
        }
    }
    else
    {
        seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf;
        if (args.samplewise)
        {
            if (estimate_args.normalization_method)
                estimate<seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>, true, true>(
                    args,
                    ibf,
                    args.path_out,
                    estimate_args);
            else
                estimate<seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>, true>(args,
                                                                                                    ibf,
                                                                                                    args.path_out,
                                                                                                    estimate_args);
        }
        else
        {
            estimate<seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>, false>(args,
                                                                                                 ibf,
                                                                                                 args.path_out,
                                                                                                 estimate_args);
        }
    }
}
