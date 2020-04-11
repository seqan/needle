#pragma once

#include <algorithm>
#include <chrono>
#include <deque>
#include <iostream>
#include <math.h>
#include <numeric>
#include <string>

#if SEQAN3_WITH_CEREAL
#include <cereal/archives/binary.hpp>
#endif // SEQAN3_WITH_CEREAL

#include <robin_hood.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/io/stream/iterator.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/std/ranges>
#include <seqan3/range/views/take_until.hpp>

#include "minimizer.h"
#include "search.h"

// specific arguments needed for constructing an IBF
struct ibf_arguments
{
    std::vector<std::filesystem::path> sequence_files;
    std::filesystem::path genome_file; // Needs to be defined when normalization_method should only use the genome sequence
    std::vector<size_t> bin_size{}; // The bin size of one IBF, can be different for different expression levels
    size_t num_hash{1}; // Number of hash functions to use, default 1
    std::filesystem::path path_out{"./"}; // Path where ibf should be stored
    std::vector<float> expression_levels{}; // 0.5,1,2,4 are default, expression levels which should be created
    std::vector<int> samples{}; // Can be used to indicate that sequence files belong to the same experiment
    bool paired = false; // If true, than experiments are seen as paired-end experiments
    // Which expression values should be ignored during calculation of the normalization_method, default is zero
    std::vector<int> cutoffs{};
    std::string normalization_method{"median"}; // Method to calculate normalized expression value
    size_t random{10}; // What percentage of sequences should be used when using normalization_method random
    bool experiment_names = false; // Flag, if names of experiment should be stored in a txt file
};

struct RandomGenerator {
	int maxi;
	RandomGenerator(int max) :
			maxi(max) {
	}

	int operator()() {
		return rand() % maxi;
	}
};

// Calculate best bin size based on number of elements maximal inserted, false positive rate and number of hash functions
// See: https://hur.st/bloomfilter/
uint64_t get_bin_size(uint64_t count, float fpr, size_t num_hash)
{
    return std::ceil((count * std::log(fpr)) / std::log(1 / std::pow(2, std::log(2))));
}

void get_sequences(std::vector<std::filesystem::path> const & sequence_files,
                   seqan3::concatenated_sequences<seqan3::dna4_vector> & sequences, unsigned min = 0, unsigned max = 1)
{
    //Loads all reads from all samples of one experiment to sequences
    //TODO: If not enough memory or too many samples in one experiment, read one file record by record
    for (unsigned i = min; i < (min + max); i++)
    {
        seqan3::sequence_file_input<my_traits> fin{sequence_files[i]};
        for (auto & [seq, id, qual]: fin)
            sequences.push_back(seq);
    }
}

void get_minimizers(arguments const & args, seqan3::concatenated_sequences<seqan3::dna4_vector> const & sequences,
                    robin_hood::unordered_node_map<uint64_t, uint64_t> & hash_table,
                    robin_hood::unordered_set<uint64_t> const & genome_set_table,
                    std::filesystem::path const & genome_file = "")
{
    // Count minimizer in sequence file
    for (auto seq : sequences)
    {
        for (auto & minHash : compute_minimizer(seq, args.k, args.window_size, args.shape, args.seed))
        {
            //TODO: Use unordered_set contains function instead of find, only works in C++20
            if ((genome_file == "") || (genome_set_table.find(minHash) != genome_set_table.end()))
                hash_table[minHash]++;
        }
    }
}

// Set arguments that ibf and minimizer use
void set_arguments(arguments const & args, ibf_arguments & ibf_args,
                   seqan3::concatenated_sequences<seqan3::dna4_vector> & genome_sequences,
                   robin_hood::unordered_set<uint64_t> & genome_set_table)
{
    if (ibf_args.paired) // If paired is true, a pair is seen as one sample
        ibf_args.samples.assign(ibf_args.sequence_files.size()/2,2);
    if (ibf_args.samples.empty()) // If no samples are given and not paired, every file is seen as one experiment
        ibf_args.samples.assign(ibf_args.sequence_files.size(),1);
    if (ibf_args.cutoffs.empty()) // If no cutoffs are given, every experiment gets a cutoff of zero
        ibf_args.cutoffs.assign(ibf_args.samples.size(),0);
    // If sum of ibf_args.samples is not equal to number of files, throw error
    else if (std::accumulate(ibf_args.samples.rbegin(), ibf_args.samples.rend(), 0) != ibf_args.sequence_files.size())
        throw std::invalid_argument{"Error. Incorrect command line input for multiple-samples."};

    // Generate genome mask
    if (ibf_args.genome_file != "")
    {
		get_sequences({ibf_args.genome_file}, genome_sequences);

        // Count minimizer in sequence file
        for (auto seq : genome_sequences)
        {
            for (auto & minHash : compute_minimizer(seq, args.k, args.window_size, args.shape, args.seed))
                genome_set_table.insert(minHash);
        }
    }

    // Sort given expression rates
    sort(ibf_args.expression_levels.begin(), ibf_args.expression_levels.end());

}

// Reads a binary file minimizer creates
void read_binary(robin_hood::unordered_node_map<uint64_t, uint64_t> & hash_table, std::filesystem::path filename)
{
    std::ifstream fin;
    uint64_t minimizer;
    uint64_t minimizer_count;
    fin.open(filename, std::ios::binary);

    while(fin.read((char*)&minimizer, sizeof(minimizer)))
    {
        fin.read((char*)&minimizer_count, sizeof(minimizer_count));
        hash_table[minimizer] = minimizer_count;
    }

    fin.close();
}

// Reads one header file minimizer creates
void read_header(arguments & args, ibf_arguments & ibf_args, std::filesystem::path filename,
                 std::vector<uint64_t> & counts, uint32_t & normalized_exp_value)
{
    std::ifstream fin;
    fin.open(filename);
    auto stream_view = seqan3::views::istreambuf(fin);
    auto stream_it = std::ranges::begin(stream_view);

    // Read first line
    std::string buffer;
    std::ranges::copy(stream_view | seqan3::views::take_until_and_consume(seqan3::is_char<' '>),
                                    std::ranges::back_inserter(buffer));
    args.seed = std::stoull(buffer);
    buffer.clear();
    std::ranges::copy(stream_view | seqan3::views::take_until_and_consume(seqan3::is_char<' '>),
                                    std::ranges::back_inserter(buffer));
    args.k = std::stoi(buffer);
    buffer.clear();
    std::ranges::copy(stream_view | seqan3::views::take_until_and_consume(seqan3::is_char<' '>),
                                    std::ranges::back_inserter(buffer));
    args.window_size = std::stoi(buffer);
    buffer.clear();
    std::ranges::copy(stream_view | seqan3::views::take_until_and_consume(seqan3::is_char<' '>),
                                    std::ranges::back_inserter(buffer));
    args.shape = std::stoull(buffer);
    buffer.clear();
    std::ranges::copy(stream_view | seqan3::views::take_until_and_consume(seqan3::is_char<' '>),
                                    std::ranges::back_inserter(buffer));
    ibf_args.normalization_method = buffer;
    buffer.clear();
    std::ranges::copy(stream_view | seqan3::views::take_until_and_consume(seqan3::is_char<'\n'>),
                                    std::ranges::back_inserter(buffer));
    normalized_exp_value = std::stoi(buffer);

    // Read second line = expression levels
    do
    {
        buffer.clear();
        std::ranges::copy(stream_view | seqan3::views::take_until_or_throw(seqan3::is_char<' '> || seqan3::is_char<'\n'>),
                                        std::ranges::back_inserter(buffer));
        ibf_args.expression_levels.push_back(std::stof(buffer));
        if (*stream_it != '\n')
            ++stream_it;
    } while (*stream_it != '\n');
    ++stream_it;

    // Read third line = counts per minimizer per expression level
    do
    {
        buffer.clear();
        std::ranges::copy(stream_view | seqan3::views::take_until_or_throw(seqan3::is_char<' '>),
                                        std::ranges::back_inserter(buffer));
        counts.push_back(std::stoull(buffer));
        if (*stream_it != '\n')
            ++stream_it;
    } while (*stream_it != '\n');
    fin.close();
}

// stores compressed and uncompressed ibfs
template <class IBFType>
void store_ibf(IBFType const & ibf,
               std::filesystem::path opath)
{
    std::ofstream os{opath, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(seqan3::interleaved_bloom_filter(ibf));
}

// Calculate normalized expression value
uint32_t normalization_method(arguments const & args, ibf_arguments const & ibf_args,
                              seqan3::concatenated_sequences<seqan3::dna4_vector> const & sequences,
                              robin_hood::unordered_node_map<uint64_t, uint64_t> & hash_table, unsigned cutoff)
{
    std::vector<uint32_t> counts;
    uint32_t mean; // the normalized expression value
    std::vector<uint32_t> medians;

    // Calculate normalized expression value by taking the median of medians of all reads of one experiment, Default
    if (ibf_args.normalization_method == "median")
    {
        for (auto seq : sequences)
        {
            for (auto minHash : compute_minimizer(seq, args.k, args.window_size, args.shape, args.seed))
            {
                counts.push_back(hash_table[minHash]);
            }
            std::nth_element(counts.begin(), counts.begin() + counts.size()/2, counts.end());
            // Do not consider expression values smaller or equal to given cutoff
            if (counts[counts.size()/2] > cutoff)
                medians.push_back(counts[counts.size()/2]);
            counts.clear();
        }
        // Take the median of the medians
        std::nth_element(medians.begin(), medians.begin() + medians.size()/2, medians.end());
        mean =  medians[medians.size()/2];
        medians.clear();
    }
    // Calculate normalized expression value by dividing the sum of all squares by the sum of all elements in hash table
    else if (ibf_args.normalization_method == "mean")
    {
        size_t sum_hash{0};
        size_t sum{0};
        // Determine mean expression of one sample
        for (auto & elem : hash_table)
        {
            if (elem.second > 0)
            {
                sum = sum + (elem.second * elem.second);
                sum_hash = sum_hash + elem.second;
            }
        }
        mean = sum/sum_hash;
    }
    // Calculate normalized expression value by taking median of medians of a given percentage of reads
    else if (ibf_args.normalization_method == "random")
    {
        // How many sequences should be looked at
        uint64_t random_n = (sequences.size()/100.0) * ibf_args.random;
        std::vector<int> randomPos(random_n);
        std::generate(randomPos.begin(), randomPos.end(), RandomGenerator(sequences.size()));
        for (auto pos : randomPos)
        {
            for (auto minHash : compute_minimizer(sequences[pos], args.k, args.window_size, args.shape, args.seed))
            {
                counts.push_back(hash_table[minHash]);
            }
            std::nth_element(counts.begin(), counts.begin() + counts.size()/2, counts.end());
            if (counts[counts.size()/2] > cutoff)
                medians.push_back(counts[counts.size()/2]);
            counts.clear();
        }
        // Take the median of the medians
        std::nth_element(medians.begin(), medians.begin() + medians.size()/2, medians.end());
        mean = medians[medians.size()/2];
        medians.clear();
    }

    return mean;
}

// Calculates statistics from header files created by minimizer
std::vector<std::tuple<std::vector<float>, std::vector<uint64_t>>> statistics(arguments & args, ibf_arguments & ibf_args,
std::vector<std::filesystem::path> const & header_files)
{
    // for every expression level a count list of all experiments is created
    std::vector<std::vector<uint64_t>> count_all{};
    std::vector<uint64_t> normalized_exp_values;
    std::vector<float> exp_levels;
    std::vector<std::tuple<std::vector<float>, std::vector<uint64_t>>> results{};

    // For function call read_header
    ibf_args.expression_levels = {};
    std::vector<uint64_t> counts{};
    uint32_t normalized_exp_value{};

    float avg;
    uint64_t minimum;
    uint64_t median;
    uint64_t maximum;
    std::vector<float> first;
    std::vector<uint64_t> second;

    for(auto & file : header_files) // Go over every minimizer file
    {
        read_header(args, ibf_args, file, counts, normalized_exp_value);
        if (count_all.size() == 0) // is true for the very first file
        {
            exp_levels = ibf_args.expression_levels;
            count_all.assign(ibf_args.expression_levels.size(), {});
        }

        for( unsigned i = 0; i < counts.size(); ++i)
            count_all[i].push_back(counts[i]);
        normalized_exp_values.push_back(normalized_exp_value);
        ibf_args.expression_levels.clear();
        counts.clear();
    }

    for( unsigned i = 0; i < exp_levels.size(); ++i)
    {
        std::nth_element(normalized_exp_values.begin(), normalized_exp_values.begin() + normalized_exp_values.size()/2,
                         normalized_exp_values.end());
        avg = (1.0 * std::accumulate(normalized_exp_values.begin(),
                                     normalized_exp_values.end(), 0))/normalized_exp_values.size();

        minimum = *std::min_element(count_all[i].begin(), count_all[i].end());
        std::nth_element(count_all[i].begin(), count_all[i].begin() + count_all[i].size()/2, count_all[i].end());
        median = count_all[i][count_all[i].size()/2];
        maximum = *std::max_element(count_all[i].begin(), count_all[i].end());

        first = {exp_levels[i], avg};
        second = {minimum, median, maximum};
        results.push_back(std::make_tuple(first, second));
    }

    return results;

}

std::vector<uint32_t> ibf(arguments const & args, ibf_arguments & ibf_args)
{
    // Declarations
    robin_hood::unordered_node_map<uint64_t, uint64_t> hash_table{}; // Storage for minimizers
    double mean; // the normalized expression value
    std::vector<uint32_t> normal_expression_values;
    seqan3::concatenated_sequences<seqan3::dna4_vector> genome_sequences; // Storage for genome sequences
    robin_hood::unordered_set<uint64_t> genome_set_table{}; // Storage for minimizers in genome sequences
    seqan3::concatenated_sequences<seqan3::dna4_vector> sequences; // Storage for sequences in experiment files

    set_arguments(args, ibf_args, genome_sequences, genome_set_table);

    //If no expression values given, add default
    if (ibf_args.expression_levels.size() == 0)
        ibf_args.expression_levels = {0.5,1,2,4};

    // If no bin size is given or not the right amount, throw error.
    if (ibf_args.bin_size.empty())
        throw std::invalid_argument{"Error. Please give a size for the IBFs in bit."};
    else if (ibf_args.bin_size.size() == 1)
        ibf_args.bin_size.assign(ibf_args.expression_levels.size(),ibf_args.bin_size[0]);
    else if (ibf_args.bin_size.size() != ibf_args.expression_levels.size())
        throw std::invalid_argument{"Error. Length of sizes for IBFs in bin_size is not equal to length of expression "
                                    "levels."};

    // Store experiment names
    if (ibf_args.experiment_names)
    {
        std::ofstream outfile;
        outfile.open(std::string{ibf_args.path_out} + "Stored_Files.txt");
        for (unsigned i = 0; i < ibf_args.samples.size(); i++)
        {
            outfile  << ibf_args.sequence_files[std::accumulate(ibf_args.samples.begin(),
                                                ibf_args.samples.begin()+i, 0)] << "\n";
        }
        outfile.close();
    }

    // Create IBFs
    std::vector<seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>> ibfs;
    for (unsigned i = 0; i < ibf_args.expression_levels.size(); i++)
        ibfs.push_back(seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>(
					   seqan3::bin_count{ibf_args.samples.size()}, seqan3::bin_size{ibf_args.bin_size[i]},
					   seqan3::hash_function_count{ibf_args.num_hash}));

    // Add minimizers to ibf
    for (unsigned i = 0; i < ibf_args.samples.size(); i++)
    {
        get_sequences(ibf_args.sequence_files, sequences, std::accumulate(ibf_args.samples.begin(),
                                                                          ibf_args.samples.begin()+i,0),
                                                                          ibf_args.samples[i]);
        get_minimizers(args, sequences, hash_table, genome_set_table, ibf_args.genome_file);

        // Calculate normalized expression value in one experiment
        // if genome file is given, calculation are based on genome sequences
        if ((ibf_args.genome_file != "") & (i==0))
        {
            mean = normalization_method(args, ibf_args, genome_sequences, hash_table, ibf_args.cutoffs[i]);
        }
        else
        {
            genome_sequences.clear();
            mean = normalization_method(args, ibf_args, sequences, hash_table, ibf_args.cutoffs[i]);
        }

        normal_expression_values.push_back(mean);
        sequences.clear();

        // Every minimizer is stored in IBF, if it occurence divided by the mean is greater or equal expression level
        for (auto & elem : hash_table)
        {
            for (unsigned j = 0; j < ibf_args.expression_levels.size(); j++)
            {
                if ((ibf_args.expression_levels[j] == 0) & (elem.second > ibf_args.cutoffs[i])) // for comparison with mantis, SBT
                    ibfs[j].emplace(elem.first,seqan3::bin_index{i});
                else if ((((double) elem.second/mean)) >= ibf_args.expression_levels[j])
                    ibfs[j].emplace(elem.first,seqan3::bin_index{i});
                else //If elem is not expressed at this level, it won't be expressed at a higher level
                    break;
            }
        }
        hash_table.clear();
    }

    // Store IBFs
    for (unsigned i = 0; i < ibf_args.expression_levels.size(); i++)
    {
        std::filesystem::path filename{ibf_args.path_out.string() + "IBF_" + std::to_string(ibf_args.expression_levels[i])};
        if (args.compressed)
        {
            seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf{ibfs[i]};
		    store_ibf(ibf, filename);
        }
        else
        {
            store_ibf(ibfs[i], filename);
        }

    }

	return normal_expression_values;

}

// Create ibf based on the minimizer and header files
std::vector<uint32_t> ibf(std::vector<std::filesystem::path> minimizer_files, std::filesystem::path header_file,
                          arguments & args, ibf_arguments & ibf_args, float fpr = 0.05)
{

    // Declarations
    robin_hood::unordered_node_map<uint64_t, uint64_t> hash_table{}; // Storage for minimizers
    seqan3::concatenated_sequences<seqan3::dna4_vector> genome_sequences; // Storage for sequences in genome
    uint32_t mean; // the normalized expression value
    std::vector<uint32_t> normal_expression_values;
    // For header file reading
    std::vector<uint64_t> counts{};
    //float normalized_exp_value{};
    std::vector<std::tuple<std::vector<float>, std::vector<uint64_t>>> statistic_results;
    // Store what user might have entered
    std::vector<float> expression_levels = ibf_args.expression_levels;
    std::string normalization_method = ibf_args.normalization_method;

    if (header_file == "") // If all header files should be considered
    {
        std::vector<std::filesystem::path> header_files{};

        for (auto & elem : minimizer_files)
            header_files.push_back(elem.replace_extension(".header"));

        statistic_results = statistics(args, ibf_args, header_files);

        for (unsigned i = 0; i < statistic_results.size(); ++i)
            counts.push_back(std::get<1>(statistic_results[i])[2]); // Get maximum count of all files

        if (expression_levels.size() == 0) // If no expression levels are given
        {
            for (unsigned i = 0; i < statistic_results.size(); ++i)
                expression_levels.push_back(std::get<0>(statistic_results[i])[0]); // Get expression levels
        }
    }
    else // Only given header file is read in
    {
        read_header(args, ibf_args, header_file, counts, mean);

        if (expression_levels.size() == 0) // If no expression levels are given
            expression_levels = ibf_args.expression_levels;
    }

    // Create IBFs
    std::vector<seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>> ibfs;
    // TODO: if expression levels given does not match the expression levels in header file, get_bin_size gives
    // incorrect results or even an error, if more expression levels are added
    for (unsigned i = 0; i < expression_levels.size(); i++)
        ibfs.push_back(seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>(
					   seqan3::bin_count{minimizer_files.size()}, seqan3::bin_size{get_bin_size(counts[i], fpr,
                                                                                                ibf_args.num_hash)},
					   seqan3::hash_function_count{ibf_args.num_hash}));

    // TODO: Take cutoffs from user?
    ibf_args.cutoffs.assign(minimizer_files.size(),0);

    // Add minimizers to ibf
    for (unsigned i = 0; i < minimizer_files.size(); i++)
    {
        read_binary(hash_table, minimizer_files[i]);

        // Get normalized expression value from header file or recalculate it when other method is asked for
        if ((normalization_method == ibf_args.normalization_method) | normalization_method == "")
        {
            read_header(args, ibf_args, minimizer_files[i].replace_extension(".header"), counts, mean);
        }
        // TODO: Add sequences to use - only genome sequences ???
        /*else
        {
            mean = normalization_method(args, ibf_args, sequences, hash_table, ibf_args.cutoffs[i]);
        }
        sequences.clear();*/

        normal_expression_values.push_back(mean);

        // Every minimizer is stored in IBF, if it occurence divided by the mean is greater or equal expression level
        for (auto & elem : hash_table)
        {
            for (unsigned j = 0; j < expression_levels.size(); j++)
            {
                if ((expression_levels[j] == 0) & (elem.second > ibf_args.cutoffs[i])) // for comparison with mantis, SBT
                    ibfs[j].emplace(elem.first,seqan3::bin_index{i});
                else if ((((double) elem.second/mean)) >= expression_levels[j])
                    ibfs[j].emplace(elem.first,seqan3::bin_index{i});
                else //If elem is not expressed at this level, it won't be expressed at a higher level
                    break;
            }
        }
        hash_table.clear();
    }

    // Store IBFs
    for (unsigned i = 0; i < expression_levels.size(); i++)
    {
        std::filesystem::path filename{ibf_args.path_out.string() + "IBF_" + std::to_string(ibf_args.expression_levels[i])};
        if (args.compressed)
        {
            seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf{ibfs[i]};
		    store_ibf(ibf, filename);
        }
        else
        {
            store_ibf(ibfs[i], filename);
        }

    }

	return normal_expression_values;
}

std::vector<uint32_t> insert(arguments const & args, ibf_arguments & ibf_args, std::filesystem::path path_in)
{
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf;
    robin_hood::unordered_node_map<uint64_t, uint64_t> hash_table{}; // Storage for minimizers
    double mean; // the normalized expression value
    std::vector<uint32_t> normal_expression_values;
    seqan3::concatenated_sequences<seqan3::dna4_vector> genome_sequences; // Storage for genome sequences
    robin_hood::unordered_set<uint64_t> genome_set_table{}; // Storage for minimizers in genome sequences
    seqan3::concatenated_sequences<seqan3::dna4_vector> sequences; // Storage for sequences in experiment files
    std::vector<size_t> bin_before; // how many bins the ibf had before

    set_arguments(args, ibf_args, genome_sequences, genome_set_table);

    // Create IBFs
    std::vector<seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>> ibfs;
    for (unsigned i = 0; i < ibf_args.expression_levels.size(); i++)
    {
        load_ibf(ibf, path_in.string() + "IBF_" + std::to_string(ibf_args.expression_levels[i]));
        bin_before.push_back(ibf.bin_count());
        ibf.increase_bin_number_to(seqan3::bin_count{ibf.bin_count() + ibf_args.samples.size()});
        ibfs.push_back(ibf);
    }

    // Add minimizers to ibf
    for (unsigned i = 0; i < ibf_args.samples.size(); i++)
    {
        get_sequences(ibf_args.sequence_files, sequences, std::accumulate(ibf_args.samples.begin(),
                                                                          ibf_args.samples.begin()+i,0),
                                                                          ibf_args.samples[i]);
        get_minimizers(args, sequences, hash_table, genome_set_table, ibf_args.genome_file);

        // Calculate normalized expression value in one experiment
        // if genome file is given, calculation are based on genome sequences
        if ((ibf_args.genome_file != "") & (i==0))
        {
            mean = normalization_method(args, ibf_args, genome_sequences, hash_table, ibf_args.cutoffs[i]);
        }
        else
        {
            genome_sequences.clear();
            mean = normalization_method(args, ibf_args, sequences, hash_table, ibf_args.cutoffs[i]);
        }

        normal_expression_values.push_back(mean);
        sequences.clear();

        // Every minimizer is stored in IBF, if it occurence divided by the mean is greater or equal expression level
        for (auto & elem : hash_table)
        {
            for (unsigned j = 0; j < ibf_args.expression_levels.size(); j++)
            {
                if ((ibf_args.expression_levels[j] == 0) & (elem.second > ibf_args.cutoffs[i])) // for comparison with mantis, SBT
                    ibfs[j].emplace(elem.first,seqan3::bin_index{bin_before[j] + i});
                else if ((((double) elem.second/mean)) >= ibf_args.expression_levels[j])
                    ibfs[j].emplace(elem.first,seqan3::bin_index{bin_before[j] + i});
                else //If elem is not expressed at this level, it won't be expressed at a higher level
                    break;
            }
        }
        hash_table.clear();
    }

    // Store IBFs
    for (unsigned i = 0; i < ibf_args.expression_levels.size(); i++)
    {
        std::filesystem::path filename{ibf_args.path_out.string() + "IBF_" + std::to_string(ibf_args.expression_levels[i])};
        store_ibf(ibfs[i], filename);
    }

    return normal_expression_values;
}

void minimizer(arguments const & args, ibf_arguments & ibf_args)
{
    // Declarations
    std::vector<uint32_t> counts;
    robin_hood::unordered_node_map<uint64_t,uint64_t> hash_table{}; // Storage for minimizers
    std::ofstream outfile;
    double mean; // the normalized expression value
    seqan3::concatenated_sequences<seqan3::dna4_vector> genome_sequences; // Storage for genome sequences
    robin_hood::unordered_set<uint64_t> genome_set_table{}; // Storage for minimizers in genome sequences
    seqan3::concatenated_sequences<seqan3::dna4_vector> sequences; // Storage for sequences in experiment files
    int seen_before{0}; // just to keep track, which sequence files have already been processed, used in stead of:
                        // std::accumulate(ibf_args.samples.begin(), ibf_args.samples.begin()+i,0)

    set_arguments(args, ibf_args, genome_sequences, genome_set_table);

    //If no expression values given, set it to zero
    if (ibf_args.expression_levels.size() == 0)
        ibf_args.expression_levels = {0};

    // Add minimizers to ibf
    for (unsigned i = 0; i < ibf_args.samples.size(); i++)
    {
        get_sequences(ibf_args.sequence_files, sequences, seen_before, ibf_args.samples[i]);
        get_minimizers(args, sequences, hash_table, genome_set_table, ibf_args.genome_file);

        // Calculate normalized expression value in one experiment
        // if genome file is given, calculation are based on genome sequences
        if ((ibf_args.genome_file != "") & (i==0))
        {
            mean = normalization_method(args, ibf_args, genome_sequences, hash_table, ibf_args.cutoffs[i]);
        }
        else
        {
            genome_sequences.clear();
            mean = normalization_method(args, ibf_args, sequences, hash_table, ibf_args.cutoffs[i]);
        }
        counts.assign(ibf_args.expression_levels.size(),0);
        for (auto & elem : hash_table)
        {
            for (unsigned j = 0; j < ibf_args.expression_levels.size(); j++)
            {
                if ((ibf_args.expression_levels[j] == 0) & (elem.second > ibf_args.cutoffs[i])) // for comparison with mantis, SBT
                    counts[j]++;
                else if ((((double) elem.second/mean)) >= ibf_args.expression_levels[j])
                    counts[j]++;
                else //If elem is not expressed at this level, it won't be expressed at a higher level
                    break;
            }
        }
        sequences.clear();

        // Write minimizer and their counts to binary
        outfile.open(std::string{ibf_args.path_out} + std::string{ibf_args.sequence_files[seen_before].stem()}
                     + ".minimizer", std::ios::binary);
        for (auto & elem : hash_table)
        {
            outfile.write((char*) &elem.first, sizeof(elem.first));
            outfile.write((char*) &elem.second, sizeof(elem.second));
        }
        outfile.close();
        hash_table.clear();

        // Write header file, containing information about the minimizer counts per expression level
        outfile.open(std::string{ibf_args.path_out} + std::string{ibf_args.sequence_files[seen_before].stem()}
                     + ".header");
        outfile << args.seed << " " << std::to_string(args.k) << " " << args.window_size << " " << args.shape << " "
                << ibf_args.normalization_method << " " << mean << "\n";
        for (unsigned k = 0; k < counts.size(); k++)
            outfile  << ibf_args.expression_levels[k] << " ";

        outfile << "\n";
        for (unsigned k = 0; k < counts.size(); k++)
            outfile  << counts[k] << " ";
        outfile << "\n";
        outfile.close();
        seen_before = seen_before + ibf_args.samples[i];
    }

}

void build_ibf(arguments & args, ibf_arguments & ibf_args, float fpr = 0.05)
{
    std::vector<std::filesystem::path> minimizer_files;
    minimizer(args, ibf_args);
    for (const auto & entry : std::filesystem::directory_iterator(ibf_args.path_out))
    {
        if (entry.path().extension() == ".minimizer")
            minimizer_files.push_back(entry.path());
    }
    ibf(minimizer_files, "", args, ibf_args, fpr);
    minimizer_files.clear();
}

void test(arguments & args, ibf_arguments & ibf_args, float fpr = 0.05, bool print = false)
{
    std::vector<std::string> methods{"median", "mean"};
    std::filesystem::path genome_file = ibf_args.genome_file;
    std::filesystem::path path_out = ibf_args.path_out;
    for(auto m: methods)
    {
        std::filesystem::create_directory(path_out/m);
        std::filesystem::create_directory(std::string(path_out/"Genome_")+m);
        auto start = std::chrono::high_resolution_clock::now();
        ibf_args.path_out = std::string(path_out/"Genome_")+m +"/";
        ibf_args.normalization_method = m;
        build_ibf(args, ibf_args, fpr);
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        if (print)
            std::cout << m << " With Genome Time taken by function: " << duration.count() << " microseconds\n";

        start = std::chrono::high_resolution_clock::now();
        ibf_args.path_out =  std::string(path_out/m) +"/";
        ibf_args.genome_file = "";
        build_ibf(args, ibf_args, fpr);
        stop = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        if (print)
            std::cout << m <<" Without Time taken by function: " << duration.count() << " microseconds\n";
    }

}
