#pragma once

#include <algorithm>
#include <deque>
#include <iostream>
#include <math.h>
#include <numeric>

#if SEQAN3_WITH_CEREAL
#include <cereal/archives/binary.hpp>
#endif // SEQAN3_WITH_CEREAL

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/std/ranges>

#include "minimizer.h"

// specific arguments needed for the search
struct search_arguments
{
    std::filesystem::path search_file;
    std::filesystem::path exp_file;
    std::filesystem::path path_in{"./"};
    float expression{1.0};

};

template <typename number_type, typename range_type>
number_type to_number(range_type && range)
{
    std::string str;
    number_type num;
    std::ranges::copy(range, std::back_inserter(str));
    auto res = std::from_chars(&str[0], &str[0] + str.size(), num);
    if (res.ec != std::errc{})
    {
        seqan3::debug_stream << "Could not cast '" << range << "' to a valid number\n";
        throw std::invalid_argument{"CAST ERROR"};
    }
    return num;
}

std::vector<uint32_t> search(arguments const & args, search_arguments const & search_args)
{
    std::vector<uint32_t> counter;
    std::vector<uint32_t> results;
    std::vector<float> expression;
    std::vector<seqan3::dna4_vector> seqs;

    seqan3::sequence_file_input<my_traits> fin{search_args.search_file};
    for (auto & [seq, id, qual] : fin)
    {
        seqs.push_back(seq);
    }

    if (search_args.exp_file != "")
    {
        throw std::invalid_argument{"Error! Number of given expression levels do not match number of sequences."};
    /*    std::ifstream file{search_args.exp_file.string()};
        if (file.is_open())
        {
            std::string line;
            while (std::getline(file, line))
            {
                auto splitted_line = line | std::view::split('\t');
                auto it = splitted_line.begin(); // move to 1rst column
                expression.push_back(to_number<double>(*std::next(it, 1)));
            }

            if (expression.size() != seqs.size())
                throw std::invalid_argument{"Error! Number of given expression levels do not match number of sequences."};
        }*/
    }
    else
    {
        expression.assign(seqs.size(),search_args.expression);
    }

    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf;
    if (args.compressed)
        seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;

    std::ifstream is{search_args.path_in.string() + "IBF_" + std::to_string(expression[0]), std::ios::binary};
    //seqan3::debug_stream << "IBF_" + std::to_string(expression[0])<< "\n";
    cereal::BinaryInputArchive iarchive{is};
    iarchive(ibf);

    uint32_t minimizer_length;
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

        minimizer_length = 0;
        for (auto & minHash : compute_minimizer(seqs[i], args.k, args.window_size, args.shape, args.seed))
        {
            std::transform (counter.begin(), counter.end(), ibf.bulk_contains(minHash).begin(), counter.begin(),
                            std::plus<int>());
            ++minimizer_length;
        }

        for(unsigned j = 0; j < counter.size(); j++)
        {
            if (counter[j] >= minimizer_length/2.0)
                results[j] = results[j] + 1;
        }
        counter.clear();
        counter.assign(ibf.bin_count(), 0);
    }

    return results;
}
