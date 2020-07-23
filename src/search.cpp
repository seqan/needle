#include <deque>
#include <iostream>
#include <math.h>
#include <numeric>
#include <string>
#include <algorithm>

#if SEQAN3_WITH_CEREAL
#include <cereal/archives/binary.hpp>
#endif // SEQAN3_WITH_CEREAL

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/std/ranges>

#include "search.h"

// Actual search
template <class IBFType>
std::vector<uint32_t> do_search(IBFType & ibf, arguments const & args, search_arguments const & search_args)
{
    std::vector<uint32_t> counter;
    std::vector<float> expression;
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
