#include <algorithm>
#include <deque>
#include <iostream>
#include <math.h>
#include <numeric>

#if SEQAN3_WITH_CEREAL
#include <cereal/archives/binary.hpp>
#endif // SEQAN3_WITH_CEREAL

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/dream_index/binning_directory.hpp>
#include <seqan3/search/dream_index/concept.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/std/ranges>

#include "minimizer3.h"

struct RandomGenerator {
	int maxi;
	RandomGenerator(int max) :
			maxi(max) {
	}

	int operator()() {
		return rand() % maxi;
	}
};

void initialize_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "Mitra Darvish";
    parser.info.short_description = "Constructs IBF on Minimizer.";
    parser.add_positional_option(args.sequence_files, "Please provide at least one sequence file.");
    parser.add_option(args.genome_file, 'g', "genom-mask", "Genom file used as a mask.");
    parser.add_option(args.bits, 'z', "size", "List of sizes in bits for IBF per expression rate.");
    parser.add_option(args.num_hash, 'n', "hash", "Number of hash functions that should be used when constructing "
                      "one IBF.");
    parser.add_option(args.path_out, 'o', "out", "Directory, where output files should be saved.");
    parser.add_option(args.expression_levels, 'e', "expression_levels", "Which expression levels should be used for "
                      "constructing the IBFs.");
    parser.add_option(args.samples, 'm', "multiple-samples", "Define which samples belong together, sum has to be equal"
                      " to number of sequence files. Default: Every sequence file is one sample from one experiment.");
    parser.add_option(args.aggregate_by, 'a', "aggregate-by", "Choose your method of aggregation: mean, median or "
                      "random. Default: median.");
    parser.add_option(args.random, 'r', "random-samples", "Choose the number of random sequences to pick from when "
                      "using aggregation method random. Default: 1000.");

    parser.add_flag(args.compressed, 'c', "compressed", "If set ibf is compressed. Default: Not compressed.");
    parser.add_option(args.k, 'k', "kmer", "Define kmer size.");
    parser.add_option(args.window_size, 'w', "window", "Define window size.");
    parser.add_option(args.shape, 'p', "shape", "Define a shape by the decimal of a bitvector, where 0 symbolizes a "
                      "position to be ignored, 1 a position considered. Default: ungapped.");
    parser.add_option(args.seed, 's', "seed", "Define seed.");
}

int main(int const argc, char const ** argv)
{

    seqan3::argument_parser parser("needle-IBF", argc, argv);
    cmd_arguments args{};
    initialize_argument_parser(parser, args);

    try
    {
        parser.parse();
    }
    catch (seqan3::parser_invalid_argument const & ext)                     // catch user errors
    {
        seqan3::debug_stream << "Error. Incorrect command line input for IBF construct." << ext.what() << "\n";
        return -1;
    }

    if (args.samples.empty()) // If no samples are given, every file is seen as on experiment
    {
        args.samples.assign(args.sequence_files.size(),1);
    }
    // If sum of args.samples is not equal to number of files
    else if (std::accumulate(args.samples.rbegin(), args.samples.rend(), 0) != args.sequence_files.size())
    {
        seqan3::debug_stream << "Error. Incorrect command line input for multiple-samples." << "\n";
        return -1;
    }

    using seqan3::get;

    // Declarations
    std::vector<uint32_t> counts;
    uint32_t file_seen{0};
    std::unordered_map<uint64_t, float> hash_table{}; // Storage for minimizers
    double mean;
    std::vector<uint32_t> medians;
    seqan3::concatenated_sequences<seqan3::dna4_vector> genome_sequences;
    seqan3::concatenated_sequences<seqan3::dna4_vector> sequences; // Storage for sequences in experiment files
    seqan3::concatenated_sequences<seqan3::dna4_vector> *sequences_ptr;

    // Sort given expression rates
    sort(args.expression_levels.begin(), args.expression_levels.end());

    //If no expression values given, add default
    if (args.expression_levels.size() == 0)
        args.expression_levels = {0.5,1,2,4};

    // If no size in bits is given or not the right amount, throw error.
    if (args.bits.empty())
    {
        seqan3::debug_stream << "Error. Please give a size for the IBFs in bit." << "\n";
        return -1;
    }
    else if (args.bits.size() == 1)
    {
        args.bits.assign(args.expression_levels.size(),args.bits[0]);
    }
    else if (args.bits.size() != args.expression_levels.size())
    {
        seqan3::debug_stream << "Error. Length of sizes for IBFs in bits is not equal to length of expression levels."
        << "\n";
        return -1;
    }

    // Generate genome mask
    std::unordered_set<uint64_t> genome_set_table{};
    if (args.genome_file != "")
    {
        seqan3::sequence_file_input<my_traits> input_file{args.genome_file};

        genome_sequences.insert(genome_sequences.end(),get<seqan3::field::SEQ>(input_file).begin(),
                                get<seqan3::field::SEQ>(input_file).end());
        // Count minimizer in fasta file
        for (auto seq : genome_sequences)
        {
            for (auto & minHash : compute_minimizer(seq, args.k, args.window_size, args.shape, args.seed))
                genome_set_table.insert(minHash);
        }
    }
    seqan3::debug_stream << "Kmers in Genome: " << genome_set_table.size() << "\n";

    // Create binning_directory
    std::vector<seqan3::binning_directory> bds;
    for (unsigned i = 0; i < args.expression_levels.size(); i++)
        bds.push_back(seqan3::binning_directory(args.samples.size(), args.bits[i], args.num_hash));

    // Add minimizers to binning_directory
    for (unsigned i = 0; i < args.samples.size(); i++)
    {
        //Loads all reads from all samples of one experiment to sequences
        //TODO: If not enough memory or too many samples in one experiment, read one file record by record
        for (unsigned ii = 0; ii < args.samples[i]; ii++)
        {
            seqan3::sequence_file_input<my_traits> input_file{args.sequence_files[file_seen]};
            sequences.insert(sequences.end(), get<seqan3::field::SEQ>(input_file).begin(),
                             get<seqan3::field::SEQ>(input_file).end());
            file_seen++;
        }
        // Count minimizer in fasta file
        for (auto seq : sequences)
        {
            for (auto & minHash : compute_minimizer(seq, args.k, args.window_size, args.shape, args.seed))
            {
                //TODO: Use unordered_set contains function instead of find, only works in C++20
                if ((args.genome_file == "") || (genome_set_table.find(minHash) != genome_set_table.end()))
                    hash_table[minHash] = hash_table[minHash] + 1.0/args.samples[i];
            }
        }

        // Calculate mean expression in one experiment
        if ((args.genome_file != "") & (i==0))
        {
            sequences.clear();
            //sequences = genome_sequences;
            sequences_ptr = &genome_sequences;
        }
        else if (i==0)
        {
            genome_sequences.clear();
            sequences_ptr = &sequences;
        }

        // Calculate mean expression by taking median of medians of all reads of one experiment, Default
        if (args.aggregate_by == "median")
        {
            for (auto seq : *sequences_ptr)
            {
                for (auto minHash : compute_minimizer(seq, args.k, args.window_size, args.shape, args.seed))
                {
                    counts.push_back(hash_table[minHash]);
                }
                std::nth_element(counts.begin(), counts.begin() + counts.size()/2, counts.end());
		// Do not consider expression values of 0
                if (counts[counts.size()/2] > 0)
                    medians.push_back(counts[counts.size()/2]);
                counts.clear();
            }

            std::nth_element(medians.begin(), medians.begin() + medians.size()/2, medians.end());
            mean =  medians[medians.size()/2];
            medians.clear();
        }
        // Calculate mean expression by dividing the sum of all squares by the sum of all elements in hash table
        else if (args.aggregate_by == "mean")
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
        // Calculate mean expression by taking median of medians of 1000 random reads of one experiment, Default
        else if (args.aggregate_by == "random")
        {
            // How many sequences should be looked at
            uint64_t random_n = (sequences.size()/100) * args.random;
            std::vector<int> randomPos(random_n);
            std::generate(randomPos.begin(), randomPos.end(), RandomGenerator(sequences_ptr->size()));
            for (auto pos : randomPos)
            {
                for (auto minHash : compute_minimizer(sequences_ptr->at(pos), args.k, args.window_size, args.shape,
                     args.seed))
                {
                    counts.push_back(hash_table[minHash]);
                }
                std::nth_element(counts.begin(), counts.begin() + counts.size()/2, counts.end());
                medians.push_back(counts[counts.size()/2]);
                counts.clear();
            }
            std::nth_element(medians.begin(), medians.begin() + medians.size()/2, medians.end());
            mean = medians[medians.size()/2];
            medians.clear();
        }
        sequences.clear();

        // Every minimizer is stored in IBF, if it occurence divided by the mean is greater or equal expression level
        for (auto & elem : hash_table)
        {
            for (unsigned j = 0; j < args.expression_levels.size(); j++)
            {
                if ((((double) elem.second/mean)) >= args.expression_levels[j])
                    bds[j].set(elem.first,i);
                else //If elem is not expressed at this level, it won't be expressed at a higher level
                    break;
            }
        }
        hash_table.clear();
    }


    // Store binning_directories
    for (unsigned i = 0; i < args.expression_levels.size(); i++)
    {
        std::filesystem::path filename{args.path_out.string() + "IBF_" + std::to_string(args.expression_levels[i])};
        std::ofstream os{filename, std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        if (args.compressed)
            oarchive(seqan3::binning_directory_compressed(bds[i]));
        else
            oarchive(seqan3::binning_directory(bds[i]));
    }

}
