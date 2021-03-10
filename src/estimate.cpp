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
std::vector<uint32_t> check_ibf(arguments const & args, IBFType & ibf, std::vector<uint32_t> & counter,
                                seqan3::dna4_vector const seq, float threshold, std::vector<uint32_t> & prev_counts,
                                std::vector<std::vector<uint32_t>> & expressions, bool last_exp, int k)
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
    for(int j = 0; j < counter.size(); j++)
    {
        if (counter[j] >= (minimiser_length * threshold))
        {
            prev_counts[j] = prev_counts[j] + counter[j];

            // Last expression level is looked at
            if (last_exp)
            {
        		results[j] = expressions[j][k];
            }
        }
        else if (prev_counts[j] >= (minimiser_length * threshold))
        {
            // Actually calculate estimation
            results[j] = expressions[j][k] - ((((minimiser_length/2.0) - counter[j])/(prev_counts[j] - (counter[j] * 1.0))) * (expressions[j][k]-expressions[j][k-1]));
        }
    }
    return results;
}

template <class IBFType>
void estimate(arguments const & args, estimate_arguments const & estimate_args, IBFType & ibf, std::filesystem::path file_out,
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
    for (auto & expression : estimate_args.expressions)
    {
        if (expression == estimate_args.expressions[estimate_args.expressions.size()-1])
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

            results = check_ibf(args, ibf, counter, seqs[i], estimate_args.threshold, prev_counts[i],
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

// Reads one header file minimiser creates
void read_header(std::vector<std::vector<uint32_t>> & expressions, std::filesystem::path filename, uint64_t iterator)
{
    std::ifstream fin;
    fin.open(filename);
    auto stream_view = seqan3::views::istreambuf(fin);
    auto stream_it = std::ranges::begin(stream_view);
    int j{0};

    // Consume first line
    std::string buffer{};
    seqan3::detail::consume(stream_view | seqan3::views::take_line_or_throw);

    // Read second line = expression levels
    do
    {
        std::ranges::copy(stream_view | seqan3::views::take_until_or_throw(seqan3::is_char<' '> || seqan3::is_char<'\n'>),
                                        std::cpp20::back_inserter(buffer));
        expressions[iterator][j] = (uint32_t)  std::stoi(buffer);
        buffer.clear();
        j++;
        if (*stream_it != '\n')
            ++stream_it;
    } while (*stream_it != '\n');
    ++stream_it;
    seqan3::detail::consume(stream_view | seqan3::views::take_line_or_throw);

    fin.close();
}

template <class IBFType>
void estimate(arguments const & args, estimate_arguments const & estimate_args, IBFType & ibf, std::filesystem::path file_out,
              std::filesystem::path search_file, std::filesystem::path path_in, std::vector<std::filesystem::path> header_files)
{
    std::vector<std::string> ids;
    std::vector<seqan3::dna4_vector> seqs;
    std::vector<uint32_t> counter;

    uint64_t minimiser_length{};

    seqan3::sequence_file_input<my_traits> fin{search_file};
    for (auto & [seq, id, qual] : fin)
    {
        ids.push_back(id);
        seqs.push_back(seq);
    }

    std::vector<std::vector<uint32_t>> prev_counts;
    bool last_exp{false};
    std::vector<std::vector<uint32_t>> expressions(header_files.size(), std::vector<uint32_t> (estimate_args.expressions.size(), 0));
    for(int i = 0; i < header_files.size(); i++)
    {
        read_header(expressions, header_files[i], i);
    }
    std::vector<uint32_t> results;
    std::vector<std::vector<uint32_t>> estimations;
    for (int j = 0; j < estimate_args.expressions.size(); j++)
    {
        if (j == estimate_args.expressions.size() - 1)
        	last_exp = true;
        load_ibf(ibf, path_in.string() + "IBF_Level_" + std::to_string(j));
        for (int i = 0; i < seqs.size(); ++i)
        {
            counter.assign(ibf.bin_count(), 0);
            if (estimations.size() <= i)
            {
                estimations.push_back(counter);
                prev_counts.push_back(counter);
            }

            results = check_ibf(args, ibf, counter, seqs[i], estimate_args.threshold, prev_counts[i],
                                expressions, last_exp, j);
            for(unsigned j = 0; j < counter.size(); j++)
                estimations[i][j] = std::max(estimations[i][j], (uint32_t) results[j]);


            counter.clear();
        }
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

void call_estimate(arguments const & args, estimate_arguments const & estimate_args, std::filesystem::path file_out,
              std::filesystem::path search_file, std::filesystem::path path_in, std::vector<std::filesystem::path> header_files = {})
{
    if (args.compressed)
    {
        seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
        if (header_files.empty())
            estimate(args, estimate_args, ibf, file_out, search_file, path_in);
        else
            estimate(args, estimate_args, ibf, file_out, search_file, path_in, header_files);
    }
    else
    {
        seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf;
        if (header_files.empty())
            estimate(args, estimate_args, ibf, file_out, search_file, path_in);
        else
            estimate(args, estimate_args, ibf, file_out, search_file, path_in, header_files);
    }
}
