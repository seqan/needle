#include <fstream>
#include <iostream>
#include <chrono>
#include <numeric>
#include <algorithm>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/std/ranges>

#include "minimizer.h"

struct cmd_arguments
{
    std::vector<std::filesystem::path> sequence_files;
    std::filesystem::path path_out{"./"};
    std::vector<int> samples{};
    bool paired{false};
    uint16_t k{20};
    uint16_t window_size{60};
    uint64_t shape;
    uint64_t seed{0x8F3F73B5CF1C9ADE};
};

void initialize_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "Mitra Darvish";
    parser.info.short_description = "Counts minimizers for a given set of experiments.";
    parser.add_positional_option(args.sequence_files, "Please provide at least one sequence file.");
    parser.add_option(args.path_out, 'o', "out", "Directory, where output files should be saved.");
    parser.add_option(args.samples, 'm', "multiple-samples", "Define which samples belong together, sum has to be equal"
                      " to number of sequence files. Default: Every sequence file is one sample from one experiment.");
    parser.add_flag(args.paired, 'p', "paired", "If experiments are paired. Default: Not paired.");
    parser.add_option(args.k, 'k', "kmer", "Define kmer size.");
    parser.add_option(args.window_size, 'w', "window", "Define window size.");
    parser.add_option(args.shape, 'p', "shape", "Define a shape by the decimal of a bitvector, where 0 symbolizes a "
                      "position to be ignored, 1 a position considered. Default: ungapped.");
    parser.add_option(args.seed, 's', "seed", "Define seed.");
}

int main(int const argc, char const ** argv)
{

    seqan3::argument_parser miniparser("needle-count", argc, argv);
    cmd_arguments args{};
    initialize_argument_parser(miniparser, args);
    seqan3::concatenated_sequences<seqan3::dna4_vector> sequences; // Storage for sequences in experiment files
    std::unordered_map<uint64_t, uint8_t> hash_table{}; // Storage for minimizers

    //cutoffs and cutoff_bounds from Mantis paper
    uint16_t cutoff{50};
    std::vector<int> cutoffs{1,3,10,20};
    std::vector<uint64_t> cutoff_bounds{314572800, 524288000, 1073741824, 3221225472};
    uint64_t count{0};
    uint64_t seen_before{0};
    std::ofstream outfile;
    //std::ifstream infile;
    uint64_t filesize{0};
    using seqan3::get;

    try
    {
        miniparser.parse();
    }
    catch (seqan3::parser_invalid_argument const & ext)                     // catch user errors
    {
        return -1;
    }

    // If no samples are given and experiments are paired, two files are seen as one experiment
    if ((args.samples.empty()) & args.paired)
    {
        args.samples.assign(args.sequence_files.size()/2,2);
    }
    else if ((args.samples.empty())) // If no samples are given, every file is seen as one experiment
    {
        args.samples.assign(args.sequence_files.size(),1);
    }
    // If sum of args.samples is not equal to number of files
    else if (std::accumulate(args.samples.rbegin(), args.samples.rend(), 0) != args.sequence_files.size())
    {
        seqan3::debug_stream << "Error. Incorrect command line input for multiple-samples." << "\n";
        return -1;
    }

    for (unsigned i = 0; i < args.samples.size(); i++)
    {
        seqan3::sequence_file_input<my_traits> fin{args.sequence_files[seen_before]};

        //sequences.push_back(*fin.begin());
        for(auto & [ seq, id, qual ] : fin)
        {
            sequences.push_back(seq);
            break;
        }

        if (sequences[0].size() < 50)
        {
            seen_before = seen_before + args.samples[i];
            seqan3::debug_stream<<"Too Short: " << i << sequences[0].size() << "\n";
            sequences.clear();
            continue;
        }
        sequences.clear();
        //Loads all reads from all samples of one experiment to sequences
        //TODO: If not enough memory or too many samples in one experiment, read one file record by record
        for (unsigned ii = 0; ii < args.samples[i]; ii++)
        {
            seqan3::sequence_file_input<my_traits> fin{args.sequence_files[seen_before + ii]};
            for(auto & [ seq, id, qual ] : fin)
            {
                sequences.push_back(seq);
            }
	        seqan3::debug_stream << "Test, Seq done\n";
            /*sequences.insert(sequences.end(), get<seqan3::field::SEQ>(input_file).begin(),
                             get<seqan3::field::SEQ>(input_file).end());*/
            filesize = filesize + std::filesystem::file_size(args.sequence_files[seen_before + ii]);
        }

    	if (filesize > 5368709120)
    	{
            seqan3::debug_stream << args.sequence_files[seen_before] << " "<< filesize << "\n";
    	    sequences.clear();
    	    seen_before = seen_before + args.samples[i];
            filesize = 0;
            count = 0;
            hash_table.clear();
            continue;
        }
        filesize = filesize + filesize; // filesize is multiplied by two because Mantis cutoff were taken for fastq, we have fasta files
        for(unsigned k = 0; k < cutoff_bounds.size(); k++)
        {
            if ((filesize) <= cutoff_bounds[k])
            {
                cutoff = cutoffs[k];
                break;
            }
        }

        filesize = 0;

        seqan3::debug_stream << "Start Minimizer\n";
        // Count minimizer in fasta file
        for (auto seq : sequences)
        {
            for (auto & minHash : compute_minimizer(seq, args.k, args.window_size, args.shape, args.seed))
            {
                hash_table[minHash] = std::min<uint8_t>(126u, hash_table[minHash]+1);
            }
        }
        seqan3::debug_stream << "End Minimizer\n";
        outfile.open(std::string{args.path_out} + std::string{args.sequence_files[seen_before].stem()} + ".minimizer", std::ios::binary);
        for (auto & minHash : hash_table)
        {
            if (minHash.second > cutoff)
            {
                outfile.write((char*) &minHash.first, sizeof(minHash.first));
                count++;
            }
        }
        outfile.close();

        outfile.open(std::string{args.path_out} + "Header_" + std::string{args.sequence_files[seen_before].stem()} + ".txt");
        seqan3::debug_stream << args.seed << "\t" << args.k << "\t" << args.window_size << "\t" << count << "\n";
        outfile  << args.seed << " " << args.k << " " << args.window_size << " " << count;
        outfile.close();
        count = 0;
        cutoff = 50;
        sequences.clear();
        hash_table.clear();
        seen_before = seen_before + args.samples[i];
        filesize = 0;

        /*//Test Binary
        infile.open(std::string{args.path_out} + std::string{args.sequence_files[seen_before-1].stem()} + ".minimizer", std::ios::binary);
        if (!infile.is_open()) {
		std::cerr << "Error in open file './in.bin'" << std::endl;
		return 1;
	}
        uint64_t num;
	    infile.read((char*)&num, sizeof(num));
        infile.close();
        seqan3::debug_stream << "Test2: " <<num << "\n";*/
    }

}
