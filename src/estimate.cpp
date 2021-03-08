#include <deque>
#include <iostream>
#include <math.h>
#include <numeric>
#include <string>
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
template <class IBFType>
std::vector<uint32_t> check_ibf(arguments const & args, IBFType & ibf, std::vector<uint32_t> & counter, seqan3::dna4_vector const seq, float threshold, std::vector<uint32_t> & prev_counts, uint64_t expression,
uint64_t prev_expression, bool last_exp)
{
    uint64_t minimiser_length = 0;
    for (auto minHash : seqan3::views::minimiser_hash(seq, args.shape, args.w_size, args.s))
    {
        auto agent = ibf.membership_agent();
        std::transform (counter.begin(), counter.end(), agent.bulk_contains(minHash).begin(), counter.begin(),
                        std::plus<int>());
        ++minimiser_length;
    }

    std::vector<uint32_t> results{};
    results.assign(ibf.bin_count(), 0);
    for(unsigned j = 0; j < counter.size(); j++)
    {
        if (counter[j] >= (minimiser_length * threshold))
        {
            prev_counts[j] = prev_counts[j] + counter[j];

            // Last expression level is looked at
            if (last_exp)
            {
        		results[j] = expression;
            }
        }
        else if (prev_counts[j] >= (minimiser_length * threshold))
        {
            // Actually calculate estimation
            results[j] = expression - ((((minimiser_length/2.0) - counter[j])/(prev_counts[j] - (counter[j] * 1.0))) * (expression-prev_expression));
        }
    }
    return results;
}

template <class IBFType>
void estimate(arguments const & args, estimate_arguments const & search_args, IBFType & ibf, std::vector<uint32_t> & expressions, std::filesystem::path file_out,
              std::filesystem::path search_file, std::filesystem::path path_in)
{
    std::vector<std::string> ids;
    std::vector<seqan3::dna4_vector> seqs;

    uint64_t minimiser_length{};

    seqan3::sequence_file_input<my_traits> fin{search_file};
    for (auto & [seq, id, qual] : fin)
    {
        ids.push_back(id);
        seqs.push_back(seq);
    }

    std::vector<std::vector<uint32_t>> prev_counts;
    uint64_t prev_expression{0};
    bool last_exp{false};

    std::vector<uint32_t> counter;
    std::vector<uint32_t> results;
    std::vector<std::vector<uint32_t>> estimations;
    for (auto & expression : expressions)
    {
        if (expression == expressions[expressions.size()-1])
        	last_exp = true;
        load_ibf(ibf, path_in.string() + "IBF_" + std::to_string(expression));
        for (int i = 0; i < seqs.size(); ++i)
        {
            counter.assign(ibf.bin_count(), 0);
            if (estimations.size() <= i)
            {
                estimations.push_back(counter);
                prev_counts.push_back(counter);
            }

            results = check_ibf(args, ibf, counter, seqs[i], search_args.threshold, prev_counts[i],
                                expression, prev_expression, last_exp);
            for(unsigned j = 0; j < counter.size(); j++)
                estimations[i][j] = std::max(estimations[i][j], (uint32_t) results[j]);


            counter.clear();
        }
        prev_expression = expression;
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

void call_estimate(arguments const & args, estimate_arguments const & search_args, std::vector<uint32_t> & expressions, std::filesystem::path file_out,
              std::filesystem::path search_file, std::filesystem::path path_in)
{
    if (args.compressed)
    {
        seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
        estimate(args, search_args, ibf, expressions, file_out, search_file, path_in);
    }
    else
    {
        seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf;
        estimate(args, search_args, ibf, expressions, file_out, search_file, path_in);
    }
}
