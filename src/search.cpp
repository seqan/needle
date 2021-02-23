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

#include "search.h"


// Check if one sequence is present in a given ibf
template <class IBFType>
std::vector<uint32_t> check_ibf(arguments const & args, IBFType & ibf, std::vector<uint32_t> & counter, seqan3::dna4_vector const seq, float threshold, std::vector<uint32_t> & prev_counts, std::vector<bool> & prev_true, uint64_t expression,
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
        if ((last_exp) & (prev_true[j]))
        {
            // results[j] = expression - ((((minimiser_length/2.0) - counter[j])/(prev_counts[j] - (counter[j] * 1.0))) * (expression-prev_expression));
    		results[j] = expression;
            prev_true[j] = false;

        }
        else if (counter[j] >= (minimiser_length * threshold))
        {
            results[j] = results[j] + 1;
            prev_counts[j] = counter[j];
        }
        else if (prev_true[j])
        {
            // seqan3::debug_stream << j << " " << minimiser_length << " " << expression << " " << prev_expression << " " << counter[j] << " " << prev_counts[j] << "\n";
            results[j] = expression - ((((minimiser_length/2.0) - counter[j])/(prev_counts[j] - (counter[j] * 1.0))) * (expression-prev_expression));
            prev_true[j] = false;
        }

    }
    return results;
}

template <class IBFType>
void estimate(arguments const & args, search_arguments const & search_args, IBFType & ibf, std::vector<uint64_t> & expressions, std::filesystem::path file_out,
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
    std::vector<std::vector<bool>> prev_true;
    std::vector<bool> prev_true1{};
    uint64_t prev_expression{0};
    bool last_exp{false};

    std::vector<uint32_t> counter;
    std::vector<uint32_t> results;
    std::vector<std::vector<uint32_t>> estimations;
    for (auto & expression : expressions)
    {
        load_ibf(ibf, path_in.string() + "IBF_" + std::to_string(expression));
        for (int i = 0; i < seqs.size(); ++i)
        {
            counter.assign(ibf.bin_count(), 0);
            if (i == 0)
            {
                prev_true1.assign(ibf.bin_count(), true);
            }
            if (estimations.size() <= i)
            {
                estimations.push_back(counter);
                prev_counts.push_back(counter);
                prev_true.push_back(prev_true1);
            }

            results = check_ibf(args, ibf, counter, seqs[i], search_args.threshold, prev_counts[i], prev_true[i], expression, prev_expression, last_exp);
            for(unsigned j = 0; j < counter.size(); j++)
                estimations[i][j] = std::max(estimations[i][j], (uint32_t) results[j]);


            counter.clear();
        }
        prev_expression = expression;
        if (expression == expressions[expressions.size()-2])
        	last_exp = true;
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

// Actual search
template <class IBFType>
std::vector<uint32_t> do_search(IBFType & ibf, arguments const & args, search_arguments const & search_args)
{
    std::vector<uint32_t> counter;
    std::vector<uint64_t> expression;
    std::vector<seqan3::dna4_vector> seqs;
    std::vector<uint32_t> results;

    seqan3::sequence_file_input<my_traits> fin{search_args.search_file};
    for (auto & [seq, id, qual] : fin)
    {
        seqs.push_back(seq);
    }

    if (search_args.exp_file != "")
    {
        std::ifstream file{search_args.exp_file.string()};
        if (file.is_open())
        {
            std::string line;
            while (std::getline(file, line))
            {
                size_t exp_start = line.find_first_not_of('\t', line.find('\t', 0));
                size_t exp_end = line.find('\n', 0);
                expression.push_back(std::stod(line.substr(exp_start, exp_end-exp_start)));
            }

            if (expression.size() != seqs.size())
                throw std::invalid_argument{"Error! Number of given expression levels do not match number of sequences."};
        }
    }
    else
    {
        expression.assign(seqs.size(),search_args.expression);
    }
    load_ibf(ibf, search_args.path_in.string() + "IBF_" + std::to_string(expression[0]));
    uint32_t minimiser_length;
    counter.resize(ibf.bin_count(), 0);
    results.resize(ibf.bin_count(), 0);
    for (unsigned i = 0; i < expression.size(); i++)
    {
        // If the expression level changes a different ibf needs to be loaded.
        // In order to keep the number of changes low, it is adviseable to perform all searches in one ibf and then
        // move on to the next expression level.
        if ((i > 0) && (expression[i] != expression[i-1]))
        {
            std::ifstream is{"IBF_" + std::to_string(expression[i]), std::ios::binary};
            //seqan3::debug_stream << "IBF_" + std::to_string(expression[i])<< "\n";
            cereal::BinaryInputArchive iarchive{is};
            iarchive(ibf);
        }

        minimiser_length = 0;
        for (auto minHash : seqan3::views::minimiser_hash(seqs[i], args.shape, args.w_size, args.s))
        {
            auto agent = ibf.membership_agent();
            std::transform (counter.begin(), counter.end(), agent.bulk_contains(minHash).begin(), counter.begin(),
                            std::plus<int>());
            ++minimiser_length;
        }

        for(unsigned j = 0; j < counter.size(); j++)
        {
            if (counter[j] >= (minimiser_length * search_args.threshold))
                results[j] = results[j] + 1;
        }
        counter.clear();
        counter.assign(ibf.bin_count(), 0);
    }

    return results;
}

void call_estimate(arguments const & args, search_arguments const & search_args, std::vector<uint64_t> & expressions, std::filesystem::path file_out,
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

std::vector<uint32_t> search(arguments const & args, search_arguments const & search_args)
{
    if (args.compressed)
    {
        seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
        return do_search(ibf, args, search_args);
    }
    else
    {
        seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf;
        return do_search(ibf, args, search_args);
    }
}
