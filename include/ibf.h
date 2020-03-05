#pragma once

#include <algorithm>
#include <deque>
#include <iostream>
#include <math.h>
#include <numeric>
#include<string>

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
                    std::unordered_map<uint64_t, uint64_t> & hash_table,
                    std::filesystem::path const & genome_file = "",
                    std::unordered_set<uint64_t> const & genome_set_table = {})
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

// Set arguments that ibf and preprocess use
void set_arguments(arguments const & args, ibf_arguments & ibf_args,
                   seqan3::concatenated_sequences<seqan3::dna4_vector> & genome_sequences,
                   std::unordered_set<uint64_t> & genome_set_table)
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

// Calculate normalized expression value
uint32_t normalization_method(arguments const & args, ibf_arguments const & ibf_args,
                              seqan3::concatenated_sequences<seqan3::dna4_vector> const & sequences,
                              std::unordered_map<uint64_t, uint64_t> & hash_table, unsigned cutoff)
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

void preprocess(arguments const & args, ibf_arguments & ibf_args)
{
    // Declarations
    std::vector<uint32_t> counts;
    std::unordered_map<uint64_t,uint64_t> hash_table{}; // Storage for minimizers
    std::ofstream outfile;
    double mean; // the normalized expression value
    seqan3::concatenated_sequences<seqan3::dna4_vector> genome_sequences; // Storage for genome sequences
    std::unordered_set<uint64_t> genome_set_table{}; // Storage for minimizers in genome sequences
    seqan3::concatenated_sequences<seqan3::dna4_vector> sequences; // Storage for sequences in experiment files

    set_arguments(args, ibf_args, genome_sequences, genome_set_table);

    //If no expression values given, set it to zero
    if (ibf_args.expression_levels.size() == 0)
        ibf_args.expression_levels = {0};

    // Add minimizers to ibf
    for (unsigned i = 0; i < ibf_args.samples.size(); i++)
    {
        get_sequences(ibf_args.sequence_files, sequences, std::accumulate(ibf_args.samples.begin(),
                                                                          ibf_args.samples.begin()+i,0),
                                                                          ibf_args.samples[i]);
        get_minimizers(args, sequences, hash_table, ibf_args.genome_file, genome_set_table);

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
        outfile.open(std::string{ibf_args.path_out} + std::string{ibf_args.sequence_files[i].stem()} + ".minimizer",
                     std::ios::binary);
        for (auto & elem : hash_table)
        {
            outfile.write((char*) &elem.first, sizeof(elem.first));
            outfile.write((char*) &elem.second, sizeof(elem.second));
        }
        outfile.close();
        hash_table.clear();

        // Write header file, containing information about the minimizer counts per expression level
        outfile.open(std::string{ibf_args.path_out} + "Header_" + std::string{ibf_args.sequence_files[i].stem()} + ".txt");
        outfile << args.seed << " " << std::to_string(args.k) << " " << args.window_size << " " << ibf_args.normalization_method << "\n";
        for (unsigned k = 0; k < counts.size(); k++)
            outfile  << ibf_args.expression_levels[k] << " ";

        outfile << "\n";
        for (unsigned k = 0; k < counts.size(); k++)
            outfile  << counts[k] << " ";
        outfile.close();
    }

}

std::vector<uint32_t> ibf(arguments const & args, ibf_arguments & ibf_args)
{
    //using seqan3::get;

    // Declarations
    std::unordered_map<uint64_t, uint64_t> hash_table{}; // Storage for minimizers
    double mean; // the normalized expression value
    std::vector<uint32_t> normal_expression_values;
    seqan3::concatenated_sequences<seqan3::dna4_vector> genome_sequences; // Storage for genome sequences
    std::unordered_set<uint64_t> genome_set_table{}; // Storage for minimizers in genome sequences
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
        get_minimizers(args, sequences, hash_table, ibf_args.genome_file, genome_set_table);

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
        std::ofstream os{filename, std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
		if (args.compressed)
		{
			seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf2{ibfs[i]};
            oarchive(seqan3::interleaved_bloom_filter(ibf2));
		}
        else
		{
            oarchive(seqan3::interleaved_bloom_filter(ibfs[i]));
		}
    }

	return normal_expression_values;

}