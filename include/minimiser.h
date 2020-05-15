#pragma once
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

#include <robin_hood.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/range/views/minimiser_hash.hpp>
#include <seqan3/search/kmer_index/shape.hpp>

//!\brief arguments used for construction of the IBF as well as the search
struct arguments
{
    bool compressed = false;
    uint8_t k{20};
    window_size w_size{60};
    seqan3::shape shape;
    seed s{0x8F3F73B5CF1C9ADE};

};

//!\brief Use dna4 instead of default dna5
struct my_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
    //TODO: Should I use a bitcompressed_vector to save memory but with the disadvantage of losing speed?
    //template <typename alph>
    //using sequence_container = seqan3::bitcompressed_vector<alph>;
};


/*! \brief Function, that calculates all minimisers for a set of sequences.
 *  \param seqs        A std::vector of sequences.
 *  \param window_size The window size to use.
 *  \param shape       The shape to use.
 *  \param seed        The seed to use.
 */
auto compute_occurrences(const std::vector<seqan3::dna4_vector> & seqs, window_size window_size, seqan3::shape shape,
                         seed seed)
{
    robin_hood::unordered_map<uint64_t, uint32_t> occurring_kmers{};
    unsigned total_count{0};

    for (auto & seq : seqs)
    {
        for (auto minHash : seqan3::views::minimiser_hash(seq, shape, window_size, seed))
        {
            ++occurring_kmers[minHash];
        }
    }

    return occurring_kmers;
}
