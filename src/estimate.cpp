// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/needle/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <deque>
#include <iostream>
#include <math.h>
#include <numeric>
#include <omp.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <algorithm>

#include <ranges>

#if SEQAN3_WITH_CEREAL
#include <cereal/archives/binary.hpp>
#endif // SEQAN3_WITH_CEREAL

#include <seqan3/alphabet/container/concatenated_sequences.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/io/sequence_file/all.hpp>

#include "estimate.h"

// Actual estimation
template <class IBFType, bool last_exp, bool normalization, typename exp_t>
void check_ibf(min_arguments const & args, IBFType const & ibf, std::vector<uint16_t> & estimations_i,
               seqan3::dna4_vector const seq, std::vector<uint32_t> & prev_counts,
               exp_t const & expressions, uint16_t const k, std::vector<double> const fprs, std::vector<int> & deleted)
{
    // Check, if one expression threshold for all or individual thresholds
    static constexpr bool multiple_expressions = std::same_as<exp_t, std::vector<std::vector<uint16_t>>>;

    // Count minimisers in ibf of current level
    std::vector<uint32_t> counter;
    counter.assign(ibf.bin_count(), 0);
    uint64_t minimiser_length = 0;
    for (auto minHash : seqan3::views::minimiser_hash(seq, args.shape, args.w_size, args.s))
    {
        auto agent = ibf.membership_agent();
        std::transform (counter.begin(), counter.end(), agent.bulk_contains(minHash).begin(), counter.begin(),
                        std::plus<int>());
        ++minimiser_length;
    }

    // Defines, where the median should be
    float minimiser_pos = minimiser_length/2.0;

    // Check every experiment by going over the number of bins in the ibf.
    for(int j = 0; j < counter.size(); j++)
    {
        if (std::find(deleted.begin(), deleted.end(), j) != deleted.end())
            continue;
        // Correction by substracting the expected number of false positives
        counter[j] = std::max((double) 0.0, (double) ((counter[j]-(minimiser_length*fprs[j]))/(1.0-fprs[j])));
        // Check, if considering previously seen minimisers and minimisers found ar current level equal to or are greater
        // than the minimiser_pow, which gives the median position.
        // If ań estimation took already place (estimations_i[j]!=0), a second estimation is not performed.
        if (((prev_counts[j] + counter[j]) >= minimiser_pos) & (estimations_i[j] == 0))
        {
            // If there was no previous level, because we are looking at the last level.
            if constexpr(last_exp)
            {
                if constexpr (multiple_expressions)
                    estimations_i[j] = expressions[k][j];
                else
                    estimations_i[j] = expressions;
            }
            else
            {
               // Actually calculate estimation, in the else case k stands for the prev_expression
               if constexpr (multiple_expressions)
                   estimations_i[j] = std::max(expressions[k][j] * 1.0, expressions[k+1][j] - ((abs(minimiser_pos - prev_counts[j])/(counter[j] * 1.0)) * (expressions[k+1][j]-expressions[k][j])));
               else
                   estimations_i[j] = std::max(expressions * 1.0, k - ((abs(minimiser_pos - prev_counts[j])/(counter[j] * 1.0)) * (k-expressions)));
            }

            // Perform normalization by dividing through the threshold of the first level. Only works, if multiple expressions were used.
            if constexpr (normalization & multiple_expressions)
                estimations_i[j] = estimations_i[j]/expressions[1][j];
        }
        else
        {
            // If not found at this level, add to previous count.
            prev_counts[j] = prev_counts[j] + counter[j];
        }
    }
}

// Reads the level file ibf creates
template<typename float_or_int>
void read_levels(std::vector<std::vector<float_or_int>> & expressions, std::filesystem::path filename)
{
    std::ifstream fin;
    fin.open(filename);
    auto stream_view = seqan3::detail::istreambuf(fin);
    auto stream_it = std::ranges::begin(stream_view);
    int j{0};
    std::vector<float_or_int> empty_vector{};

    std::string buffer{};

    // Read line = expression levels
    do
    {
        if (j == expressions.size())
            expressions.push_back(empty_vector);
        std::ranges::copy(stream_view | seqan3::detail::take_until_or_throw(seqan3::is_char<' '>),
                                        std::back_inserter(buffer));
        if constexpr(std::same_as<uint16_t, float_or_int>)
            expressions[j].push_back((uint16_t)  std::stoi(buffer));
        else
            expressions[j].push_back((double)  std::stod(buffer));
        buffer.clear();
        if(*stream_it != '/')
            ++stream_it;

        if (*stream_it == '\n')
        {
            ++stream_it;
            j++;
        }
    } while (*stream_it != '/');
    ++stream_it;

    fin.close();
}

/*! \brief Function to estimate expression value.
*  \param args        The arguments.
*  \param ibf         The ibf determing what kind ibf is used (compressed or uncompressed).
*  \param file_out    The output file.
*  \param estimate_args  The estimate arguments.
*/
template <class IBFType, bool samplewise, bool normalization_method = false>
void estimate(estimate_ibf_arguments & args, IBFType & ibf, std::filesystem::path file_out,
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

    seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::id, seqan3::field::seq>> fin{estimate_args.search_file};
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
        load_ibf(ibf, estimate_args.path_in.string() + "IBF_Level_" + std::to_string(args.number_expression_thresholds-1));
    else
        load_ibf(ibf, estimate_args.path_in.string() + "IBF_" + std::to_string(args.expression_thresholds[args.expression_thresholds.size()-1]));
    counter.assign(ibf.bin_count(), 0);
    counter_est.assign(ibf.bin_count(), 0);

    for (int i = 0; i < seqs.size(); ++i)
    {
        estimations.push_back(counter_est);
        prev_counts.push_back(counter);
    }
    counter_est.clear();
    counter.clear();

    // Go over the sequences
    #pragma omp parallel for
    for (int i = 0; i < seqs.size(); ++i)
    {
        if constexpr (samplewise & normalization_method)
            check_ibf<IBFType, true, true>(args, ibf, estimations[i], seqs[i], prev_counts[i],
                                           expressions,args.number_expression_thresholds - 1,
                                           fprs[args.number_expression_thresholds - 1], deleted);
        else if constexpr (samplewise)
            check_ibf<IBFType, true, false>(args, ibf, estimations[i], seqs[i], prev_counts[i],
                                            expressions, args.number_expression_thresholds - 1,
                                            fprs[args.number_expression_thresholds - 1], deleted);
        else
            check_ibf<IBFType, true, false>(args, ibf, estimations[i], seqs[i], prev_counts[i],
                                            args.expression_thresholds[args.expression_thresholds.size() - 1], prev_expression,
                                            fprs[args.expression_thresholds.size() - 1], deleted);
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
        for (int i = 0; i < seqs.size(); ++i)
        {
            if constexpr (samplewise & normalization_method)
                check_ibf<IBFType, false, true>(args, ibf, estimations[i], seqs[i], prev_counts[i],
                                      expressions, j, fprs[j], deleted);
            else if constexpr (samplewise)
                check_ibf<IBFType, false, false>(args, ibf, estimations[i], seqs[i], prev_counts[i],
                                          expressions, j, fprs[j], deleted);
            else
                check_ibf<IBFType, false, false>(args, ibf, estimations[i], seqs[i], prev_counts[i],
                                          args.expression_thresholds[j], prev_expression, fprs[j], deleted);
        }

        if (!samplewise)
            prev_expression = args.expression_thresholds[j];
    }

    // Write output file.
    std::ofstream outfile;
    outfile.open(std::string{file_out});
    for (int i = 0; i <  seqs.size(); ++i)
    {
        outfile << ids[i] << "\t";
        for (int j = 0; j < ibf.bin_count(); ++j)
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
                estimate<seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed>, true, true>(args, ibf, args.path_out, estimate_args);
            else
                estimate<seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed>, true>(args, ibf, args.path_out, estimate_args);
        }
        else
        {
            estimate<seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed>, false>(args, ibf, args.path_out, estimate_args);
        }
    }
    else
    {
        seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf;
        if (args.samplewise)
        {
            if (estimate_args.normalization_method)
                estimate<seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>, true, true>(args, ibf, args.path_out, estimate_args);
            else
                estimate<seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>, true>(args, ibf, args.path_out, estimate_args);
        }
        else
        {
            estimate<seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>, false>(args, ibf, args.path_out, estimate_args);
        }
    }
}
