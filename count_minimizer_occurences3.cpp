#include <iostream>
#include <chrono>
#include <numeric>
#include <algorithm>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>
//#include <seqan3/io/stream/debug_stream.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/std/ranges>

#include "minimizer3.h"

uint8_t num_bins;

uint64_t constexpr num_possible_kmers{274877906944}; // 4**k
uint8_t constexpr min_num_minimizer_per_read{3}; // per read

uint8_t k;
size_t window_size;
uint64_t seed;

using namespace seqan3;
// Change default traits from sequence_file
struct my_traits : sequence_file_input_default_traits_dna
{
    using sequence_alphabet = dna4;               // instead of dna5
};

struct cmd_arguments
{
    std::filesystem::path file_path_fasta{};
    std::filesystem::path file_path_exon;
    std::filesystem::path file_path_out{"out.txt"};
    uint8_t k{20};
    uint16_t window_size{60};
    uint64_t shape;
    uint64_t seed{0x8F3F73B5CF1C9ADE};
};

void initialize_argument_parser(argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "Mitra Darvish";
    parser.info.short_description = "Calculates Minimizer.";
    parser.add_positional_option(args.file_path_fasta, "Please provide a fasta file.");
    parser.add_option(args.file_path_exon, 'e', "exon", "Fasta file with exon information.");
    parser.add_option(args.file_path_out, 'o', "out", "File to which the output is saved to.");
    parser.add_option(args.k, 'k', "kmer", "Define kmer size.");
    parser.add_option(args.window_size, 'w', "window", "Define window size.");
    parser.add_option(args.shape, 'p', "shape", "Define a shape by the decimal of a bitvector, where 0 symbolizes a "
                      "position to be ignored, 1 a position considered. Default: ungapped.");
    parser.add_option(args.seed, 's', "seed", "Define seed.");
}

/*
std::vector<uint64_t> compute_minimizer(dna4_vector & seq)
{
//    size_t const window_size = ((length(seq) - k) / LeastNumMinimizers) + k; // ensure at least 3 minimizer
    if (k > window_size)
        throw std::logic_error("Cannot divide the sequence into the specified window size."); // return std::array<uint64_t, 3> {};

    Minimizer minimizer;
    minimizer.resize(k, window_size, seed);

    return minimizer.getMinimizer(seq);
}


auto compute_occurrences(std::vector<dna4_vector> & seqs)
{
    std::unordered_map<uint64_t, uint32_t> occurring_kmers{};

    for (auto & seq : seqs)
    {
        for (auto & minHash : compute_minimizer(seq))
        {
            ++occurring_kmers[minHash];
        }
    }

    return occurring_kmers;
}
*/

int main(int const argc, char const ** argv)
{

    argument_parser miniparser("Minimizer", argc, argv);
    cmd_arguments args{};
    initialize_argument_parser(miniparser, args);

    try
    {
        miniparser.parse();
    }
    catch (parser_invalid_argument const & ext)                     // catch user errors
    {
        //debug_stream << "Error. Incorrect command line input for minimizer." << ext.what() << "\n";
        return -1;
    }

    k = args.k;
    window_size = args.window_size;
    seed = args.seed;
    char const sep{'\t'};
    std::ofstream outfile;
    outfile.open(args.file_path_out);

    std::vector<std::string> ids;
    std::vector<dna4_vector> seqs;

    std::cerr << args.file_path_fasta << std::endl;
    std::cerr << "Reading in sequences." << std::endl;
    auto start = std::chrono::steady_clock::now();
    sequence_file_input<my_traits> input_file{args.file_path_fasta};
    for (auto & rec : input_file)
    {
        ids.push_back(get<field::ID>(rec));
        seqs.push_back(get<field::SEQ>(rec));
    }
    // TODO. one can cluster at a per read basis! without reading every read into memory
    auto end = std::chrono::steady_clock::now();
    std::cerr << "Time: " << std::chrono::duration_cast<std::chrono::seconds>(end-start).count() << "s" << std::endl << std::endl;

    //window_size = ((length(seqs[0]) - k) / min_num_minimizer_per_read) + k; // ensure at least 3 minimizer
    std::cerr << "Window size used for this set: " << args.window_size << std::endl;

    std::cerr << "Computing minimizers in all sequences." << std::endl;

    start = std::chrono::steady_clock::now();
    auto hash_table = compute_occurrences(seqs, k, window_size, args.shape, seed);
    end = std::chrono::steady_clock::now();
    std::cerr << "Time: " << std::chrono::duration_cast<std::chrono::seconds>(end-start).count() << "s" << std::endl << std::endl;

    // Compute and print the median occurence for each read
    std::vector<uint32_t> counts;
    size_t seq_idx{0};
    for (auto & seq : seqs)
    {
        for (auto & minHash : compute_minimizer(seq, k, window_size, args.shape, seed))
        {
            counts.push_back(hash_table[minHash]);
            outfile << ids[seq_idx] << sep << minHash << "\n";
        }
        std::nth_element(counts.begin(), counts.begin() + counts.size()/2, counts.end());
        if ( 1 < counts[counts.size()/2] )
            debug_stream << ids[seq_idx] << sep << counts[counts.size()/2] << "\n";
        counts.clear();
        seq_idx++;
    }

    if(args.file_path_exon != "") // exon file given
    {
        // assume a exon fasta file was given
        std::vector<std::string> exon_ids;
        std::vector<dna4_vector> exon_seqs;

        sequence_file_input<my_traits> gene_file{args.file_path_exon};
        for (auto & rec : gene_file)
        {
          exon_ids.push_back(get<field::ID>(rec));
          exon_seqs.push_back(get<field::SEQ>(rec));
        }

        for (size_t idx = 0; idx < exon_ids.size(); ++idx)
        {
            //outfile << exon_ids[idx] << "\t" << idx << "\t[";
            //debug_stream << exon_ids[idx] << "\n";
            for (auto hash : compute_minimizer(exon_seqs[idx], k, window_size, args.shape, seed))
            {
            //    outfile << hash_table[hash] << ",";
                counts.push_back(hash_table[hash]);
            //    debug_stream << hash << " ";
            }
            //outfile << "]\n";
            std::nth_element(counts.begin(), counts.begin() + counts.size()/2, counts.end());
            //std::sort(counts.begin(), counts.end());
            //debug_stream << counts << "\n";
            //if (counts[counts.size()/2] > 1)
            debug_stream << exon_ids[idx] << sep << counts[counts.size()/2] << sep << idx << "\n";

            counts.clear();
        }
    }
    outfile.close();
}
