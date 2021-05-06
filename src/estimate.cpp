#include <deque>
#include <iostream>
#include <math.h>
#include <numeric>
#include <stdlib.h>
#include <string>
#include <vector>
#include <algorithm> //reorded because of this error:https://github.com/Homebrew/homebrew-core/issues/44579


#if SEQAN3_WITH_CEREAL
#include <cereal/archives/binary.hpp>
#endif // SEQAN3_WITH_CEREAL

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/std/ranges>
#include <seqan3/range/container/concatenated_sequences.hpp>

#include "estimate.h"


// Check if one sequence is present in a given ibf
template <class IBFType, bool last_exp>
void check_ibf(arguments const & args, IBFType const & ibf, std::vector<uint16_t> & estimations_i, seqan3::dna4_vector const seq,
               std::vector<uint32_t> & prev_counts, uint16_t expression, uint16_t prev_expression)
{
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

    float minimiser_pos = minimiser_length/2.0;

    for(unsigned j = 0; j < counter.size(); j++)
    {
        if ((prev_counts[j] + counter[j]) >= minimiser_pos)
        {
            if constexpr(last_exp)
            {
                estimations_i[j] = expression;
            }
            else
            {
                // Actually calculate estimation
                estimations_i[j] = expression + ((abs(minimiser_pos - counter[j])/abs((prev_counts[j]*1.0) - counter[j])) * (prev_expression-expression));
                // Make sure, every transcript is only estimated once
                prev_counts[j] = 0;
            }
        }
        else
        {
            prev_counts[j] = prev_counts[j] + counter[j];
        }
    }

}

template <class IBFType, bool last_exp>
void check_ibf(arguments const & args, IBFType const & ibf, std::vector<uint16_t> & estimations_i,
               seqan3::dna4_vector const seq, std::vector<uint32_t> & prev_counts,
               std::vector<std::vector<uint16_t>> & expressions, int k)
{
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

    float minimiser_pos = minimiser_length/2.0;

    for(int j = 0; j < counter.size(); j++)
    {
        if ((prev_counts[j] + counter[j]) >= minimiser_pos)
        {
            // If there was nothing previous
            if constexpr(last_exp)
            {
                estimations_i[j] = expressions[k][j];
            }
            else
            {
                // Actually calculate estimation
                estimations_i[j] = expressions[k][j] + ((abs(minimiser_pos - counter[j])/abs((prev_counts[j]*1.0) - counter[j])) * (expressions[k+1][j]-expressions[k][j]));
                // Make sure, every transcript is only estimated once
                prev_counts[j] = 0;
            }
        }
        else
        {
            prev_counts[j] = prev_counts[j] + counter[j];
        }
    }

}

// Reads the level file ibf creates
void read_levels(std::vector<std::vector<uint16_t>> & expressions, std::filesystem::path filename)
{
    std::ifstream fin;
    fin.open(filename);
    auto stream_view = seqan3::views::istreambuf(fin);
    auto stream_it = std::ranges::begin(stream_view);
    int j{0};
    std::vector<uint16_t> empty_vector{};

    std::string buffer{};

    // Read line = expression levels
    do
    {
        if (j == expressions.size())
            expressions.push_back(empty_vector);
        std::ranges::copy(stream_view | seqan3::views::take_until_or_throw(seqan3::is_char<' '>),
                                        std::cpp20::back_inserter(buffer));
        expressions[j].push_back((uint16_t)  std::stoi(buffer));
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
*  \param estimate_args The search arguments.
*  \param ibf         The ibf determing what kind ibf is used (compressed or uncompressed).
*  \param file_out    The file where results should be stored to.
*  \param search_file The sequence file with the sequences which expression value should be estimated.
*  \param path_in     The directory where the ibfs can be found.
*  \param level_file Path to the header files, where expression levels can be found.
*/
template <class IBFType, bool samplewise>
void estimate(arguments const & args, estimate_arguments & estimate_args, IBFType & ibf, std::filesystem::path file_out,
              std::filesystem::path search_file, std::filesystem::path path_in, std::filesystem::path level_file = "")
{
    std::vector<std::string> ids;
    std::vector<seqan3::dna4_vector> seqs;
    std::vector<uint32_t> counter;
    std::vector<uint16_t> counter_est;
    std::vector<std::vector<uint32_t>> prev_counts;
    uint64_t prev_expression;
    std::vector<std::vector<uint16_t>> estimations;
    std::vector<std::vector<uint16_t>> expressions;

    seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::id, seqan3::field::seq>> fin{search_file};
    for (auto & [id, seq] : fin)
    {
        ids.push_back(id);
        seqs.push_back(seq);
    }

    if constexpr (samplewise)
        read_levels(expressions, level_file);
    else
        prev_expression = 0;

    // Make sure expression levels are sorted.
    sort(estimate_args.expressions.begin(), estimate_args.expressions.end());

    // Initialse last expression
    if constexpr (samplewise)
        load_ibf(ibf, path_in.string() + "IBF_Level_" + std::to_string(estimate_args.expressions.size()-1));
    else
        load_ibf(ibf, path_in.string() + "IBF_" + std::to_string(estimate_args.expressions[estimate_args.expressions.size()-1]));
    counter.assign(ibf.bin_count(), 0);
    counter_est.assign(ibf.bin_count(), 0);

    for (int i = 0; i < seqs.size(); ++i)
    {
        estimations.push_back(counter_est);
        prev_counts.push_back(counter);
    }
    counter_est.clear();
    counter.clear();
    for (int i = 0; i < seqs.size(); ++i)
    {
        if constexpr (samplewise)
            check_ibf<IBFType, true>(args, ibf, estimations[i], seqs[i], prev_counts[i],
                                     expressions, estimate_args.expressions.size() - 1);
        else
            check_ibf<IBFType, true>(args, ibf, estimations[i], seqs[i], prev_counts[i],
                                     estimate_args.expressions[estimate_args.expressions.size() - 1], prev_expression);
    }

    if constexpr (!samplewise)
        prev_expression = estimate_args.expressions[estimate_args.expressions.size() - 1];

    for (int j = estimate_args.expressions.size() - 2; j >= 0; j--)
    {
        if constexpr (samplewise)
            load_ibf(ibf, path_in.string() + "IBF_Level_" + std::to_string(j));
        else
            load_ibf(ibf, path_in.string() + "IBF_" + std::to_string(estimate_args.expressions[j]));
        for (int i = 0; i < seqs.size(); ++i)
        {
            if constexpr (samplewise)
                check_ibf<IBFType, false>(args, ibf, estimations[i], seqs[i], prev_counts[i],
                                          expressions, j);
            else
                check_ibf<IBFType, false>(args, ibf, estimations[i], seqs[i], prev_counts[i],
                                          estimate_args.expressions[j], prev_expression);
        }

        if (!samplewise)
            prev_expression = estimate_args.expressions[j];
    }

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

void call_estimate(arguments const & args, estimate_arguments & estimate_args, std::filesystem::path file_out,
                   std::filesystem::path search_file, std::filesystem::path path_in, std::filesystem::path level_file)
{
    if (args.compressed)
    {
        seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
        if (level_file == "")
            estimate<seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed>, false>(args, estimate_args, ibf, file_out, search_file, path_in);
        else
            estimate<seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed>, true>(args, estimate_args, ibf, file_out, search_file, path_in, level_file);
    }
    else
    {
        seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf;
        if (level_file == "")
            estimate<seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>, false>(args, estimate_args, ibf, file_out, search_file, path_in);
        else
            estimate<seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>, true>(args, estimate_args, ibf, file_out, search_file, path_in, level_file);
    }
}
