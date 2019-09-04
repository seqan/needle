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

using namespace seqan3;

struct cmd_arguments
{
    std::vector<std::filesystem::path> fasta_files;
    //std::filesystem::path gene_file;
    std::filesystem::path genome_file;
    std::vector<size_t> bits{};
    size_t num_hash{1};
    size_t random{10};
    std::filesystem::path path_out{"./"};
    std::vector<float> expression{}; // 0.5,1,2,4
    std::vector<int> samples{};
    std::string aggregate_by{"median"};
    bool compressed = false;
    uint8_t k{20};
    uint16_t window_size{60};
    uint64_t shape;
    uint64_t seed{0x8F3F73B5CF1C9ADE};
};

// Change default traits from sequence_file
struct my_traits : sequence_file_input_default_traits_dna
{
    using sequence_alphabet = dna4;               // instead of dna5
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

void initialize_argument_parser(argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "Mitra Darvish";
    parser.info.short_description = "Constructs IBF on Minimizer.";
    parser.add_positional_option(args.fasta_files, "Please provide at least one fasta file. Every file is threated as a "
                                 "result of one RNA sequence experiment.", input_file_validator{{"fa","fasta"}});
    parser.add_option(args.bits, 'l', "size", "List of sizes in bits for IBF per expression rate.");
    parser.add_option(args.num_hash, 'n', "hash", "Number of hash functions used.");
    parser.add_option(args.path_out, 'o', "out", "Directory where output files should be saved.");
    parser.add_option(args.expression, 'e', "expression_list", "Which expression levels should be used for constructing "
                      "an IBF.");
    parser.add_option(args.k, 'k', "kmer", "Define kmer size.");
    parser.add_option(args.window_size, 'w', "window", "Define window size.");
    parser.add_option(args.shape, 'p', "shape", "Define a shape by the decimal of a bitvector, where 0 symbolizes a "
                      "position to be ignored, 1 a position considered. Default: ungapped.");
    parser.add_option(args.seed, 's', "seed", "Define seed.");
    //parser.add_option(args.gene_file, 'g', "gene", "Gene files");
    parser.add_option(args.genome_file, 'g', "genom", "Genom file used as a mask.");
    parser.add_option(args.samples, 'm', "multiple-samples", "Define which samples belong together, sum has to be equal "
                      "to number of fasta files. Default: Every fasta file is one sample from one experiment.");
    parser.add_option(args.aggregate_by, 'a', "aggregate-by", "Choose your method of aggregation: mean, median or random."
                      " Default: median.");
    parser.add_option(args.random, 'r', "random-samples", "Choose the number of random sequences to pick from when "
                      "using method random. Default: 1000.");
    parser.add_flag(args.compressed, 'c', "compressed", "If set imbf is compressed. Default: Not compressed.");
}

int main(int const argc, char const ** argv)
{

    argument_parser parser("IMBF", argc, argv);
    cmd_arguments args{};
    initialize_argument_parser(parser, args);

    try
    {
        parser.parse();
    }
    catch (parser_invalid_argument const & ext)                     // catch user errors
    {
        debug_stream << "Error. Incorrect command line input for IBF construct." << ext.what() << "\n";
        return -1;
    }

    if (args.samples.empty()) // If no samples are given, every file is seen as on experiment
    {
        args.samples.assign(args.fasta_files.size(),1);
    }
    // If sum of args.samples is not equal to number of files
    else if (std::accumulate(args.samples.rbegin(), args.samples.rend(), 0) != args.fasta_files.size())
    {
        debug_stream << "Error. Incorrect command line input for multiple-samples." << "\n";
        return -1;
    }

    // Declarations
    std::vector<uint32_t> counts;
    uint32_t file_seen{0};
    std::unordered_map<uint64_t, float> hash_table{}; // Storage for minimizers
    size_t mean;
    std::vector<uint32_t> medians;
    concatenated_sequences<dna4_vector> genome_sequences;
    concatenated_sequences<dna4_vector> sequences; // Storage for sequences in experiment files
    //TODO: Delete, just for analysis
    std::vector<uint32_t> means_list;
    size_t seq_minimizer = 0;
    std::vector<uint64_t> elem_saved_exp = {0,0,0,0,0,0};

    // Sort given expression rates
    sort(args.expression.begin(), args.expression.end());

    //If no expression values given, add default
    if (args.expression.size() == 0)
        args.expression = {0.5,1,2,4};

    // If no size in bits is given or not the right amount, throw error.
    if (args.bits.empty())
    {
        debug_stream << "Error. Please give a size for the IBFs in bit." << "\n";
        return -1;
    }
    else if (args.bits.size() == 1)
    {
        args.bits.assign(args.expression.size(),args.bits[0]);
    }
    else if (args.bits.size() != args.expression.size())
    {
        debug_stream << "Error. Length of sizes for IBFs in bits is not equal to length of expression levels. " << "\n";
        return -1;
    }

    // Generate genome mask
    std::unordered_set<uint64_t> genome_set_table{};
    if (args.genome_file != "")
    {
        sequence_file_input<my_traits> input_file{args.genome_file};
        genome_sequences.insert(genome_sequences.end(), get<field::SEQ>(input_file).begin(), get<field::SEQ>(input_file).end());
        // Count minimizer in fasta file
        for (auto seq : genome_sequences)
        {
            for (auto & minHash : compute_minimizer(seq, args.k, args.window_size, args.shape, args.seed))
                genome_set_table.insert(minHash);
        }
    }
    debug_stream << "Kmers in Genome: " << genome_set_table.size() << "\n";

    // Create binning_directory
    std::vector<binning_directory> bds;
    for (unsigned i = 0; i < args.expression.size(); i++)
        bds.push_back(binning_directory(args.samples.size(), args.bits[i], args.num_hash));

    // Add minimizers to binning_directory
    for (unsigned i = 0; i < args.samples.size(); i++)
    {
        auto start = std::chrono::steady_clock::now();
        //Loads all reads from all samples of one experiment to sequences
        //TODO: If not enough memory or too many samples in one experiment, read one file record by record
        for (unsigned ii = 0; ii < args.samples[i]; ii++)
        {
            sequence_file_input<my_traits> input_file{args.fasta_files[file_seen]};
            sequences.insert(sequences.end(), get<field::SEQ>(input_file).begin(), get<field::SEQ>(input_file).end());
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

        auto end = std::chrono::steady_clock::now();
        std::cerr << "Time, Minimizer created: " << std::chrono::duration_cast<std::chrono::seconds>(end-start).count() << "s" << std::endl << std::endl;

        //TODO: Delete, just for analysis
        for (auto & elem : hash_table)
        {
            if (elem.second > 0)
            {
                seq_minimizer++;
            }
        }

        debug_stream << "Seq minimizer: " << seq_minimizer << "\n";
        seq_minimizer = 0;

        // Calculate mean expression in one experiment
        start = std::chrono::steady_clock::now();
        if (args.genome_file != "")
        {
            sequences.clear();
            sequences = genome_sequences;
        }
        else
        {
            genome_sequences.clear();
        }

        // Calculate mean expression by taking median of medians of all reads of one experiment, Default
        if (args.aggregate_by == "median")
        {
            for (auto seq : sequences)
            {
                for (auto minHash : compute_minimizer(seq, args.k, args.window_size, args.shape, args.seed))
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
            uint64_t random_n = 1000; //{(sequences.size()/100) * args.random};
            std::vector<int> randomPos(random_n);
            std::generate(randomPos.begin(), randomPos.end(), RandomGenerator(sequences.size()));
            for (auto pos : randomPos)
            {
                for (auto minHash : compute_minimizer(sequences[pos], args.k, args.window_size, args.shape, args.seed))
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
        end = std::chrono::steady_clock::now();
        std::cerr << "Time, Mean expression calculated: " << std::chrono::duration_cast<std::chrono::seconds>(end-start).count() << "s" << std::endl << std::endl;
        sequences.clear();
        means_list.push_back(end-start);

        // Every minimizer is stored in IBF, if it occurence divided by the mean is greater or equal expression level
        for (auto & elem : hash_table)
        {
            for (unsigned j = 0; j < args.expression.size(); j++)
            {
                if ((((double) elem.second/mean) + 0.05) >= args.expression[j])
                {
                    if ( bds[j].get(elem.first)[i] == 0)
                        elem_saved_exp[j]++;
                    bds[j].set(elem.first,i);
                }
                else //If elem is not expressed at this level, it won't be expressed at a higher level
                    break;
            }
        }


        hash_table.clear();
    }

    // Store binning_directories
    for (unsigned i = 0; i < args.expression.size(); i++)
    {
        std::filesystem::path filename{args.path_out.string() + "IMBF_" + std::to_string(args.expression[i])};
        std::ofstream os{filename, std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        if (args.compressed)
            oarchive(binning_directory_compressed(bds[i]));
        else
            oarchive(binning_directory(bds[i]));
    }

    // Just for testing, TODO: delete later
    double sum = std::accumulate(means_list.begin(), means_list.end(), 0.0);
    double meani = sum / means_list.size();
    debug_stream << "Mean: "<< meani << "\n";
    double sq_sum = std::inner_product(means_list.begin(),means_list.end(), means_list.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / means_list.size() - meani * meani);
    debug_stream << "Std: "<< stdev<< "\n";
    debug_stream << "Seq minimizer: " << seq_minimizer << "\n";
    debug_stream << "Elements per expression rate: "<<elem_saved_exp << "\n";
    debug_stream << "Built IBF\n";
    debug_stream << args.fasta_files.size() << "\n";
    debug_stream << args.expression.size() << "\n";
    debug_stream << args.expression << "\n";
/*
    std::vector<dna4_vector> genes;
    sequence_file_input input_file{args.genome_file};
    std::vector<int> counter_help;
    for (unsigned i = 0; i < args.samples.size(); i++)
        counter_help.push_back(0);
    std::vector<std::vector<int>> count_genes;
    for (unsigned i = 0; i < args.expression.size(); i++)
        count_genes.push_back(counter_help);

    for (auto & rec : input_file)
        genes.push_back(get<field::SEQ>(rec));

    for (auto & seq : genes)
    {
        std::vector<std::vector<int>> counter;
        int s{0};
        for (unsigned i = 0; i < args.expression.size(); i++)
            counter.push_back(counter_help);
        for (auto & minHash : compute_minimizer(seq, args.k, args.window_size, args.seed))
        {
            for (unsigned i = 0; i < args.expression.size(); i++)
            {
                for (unsigned j = 0; j < args.samples.size(); j++)
                {
                        if (bds[i].get(minHash)[j])
                            counter[i][j]++;
                }
            }
            s++;
        }

        for (unsigned i = 0; i < args.expression.size(); i++)
        {
            for (unsigned j = 0; j < args.samples.size(); j++)
            {
                if (counter[i][j] > s/2)
                    count_genes[i][j]++;
                debug_stream << i<< j<<" "<< counter[i][j]<<" "<< s << "\n";
            }
        }
        counter.clear();
    }
    debug_stream << count_genes << "\n";*/

}
