#pragma once
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

#include <robin_hood.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/range/views/minimiser_hash.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/search/kmer_index/shape.hpp>

inline constexpr static uint64_t adjust_seed(uint8_t const kmer_size, uint64_t const seed = 0x8F3F73B5CF1C9ADEULL) noexcept
{
    return seed >> (64u - 2u * kmer_size);
}

//!\brief arguments used for construction of the IBF as well as the search
struct arguments
{
    bool compressed = false;
    uint8_t k{20};
    std::filesystem::path path_out{"./"};
    seqan3::seed s{0x8F3F73B5CF1C9ADEULL};
    seqan3::shape shape = seqan3::ungapped{k};
    uint8_t threads{1};
    seqan3::window_size w_size{60};
};

//!\brief Use dna4 instead of default dna5
struct my_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
    //TODO: Should I use a bitcompressed_vector to save memory but with the disadvantage of losing speed?
    //template <typename alph>
    //using sequence_container = seqan3::bitcompressed_vector<alph>;
};

/*! \brief Function, loading compressed and uncompressed ibfs
 *  \param ibf   ibf to load
 *  \param ipath Path, where the ibf can be found.
 */
template <class IBFType>
void load_ibf(IBFType & ibf, std::filesystem::path ipath)
{
    std::ifstream is{ipath, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    iarchive(ibf);
}

/*! \brief Function, which stored compressed and uncompressed ibfs
 *  \param ibf   The IBF to store.
 *  \param opath Path, where the IBF can be found.
 */
template <class IBFType>
void store_ibf(IBFType const & ibf,
               std::filesystem::path opath)
{
    std::ofstream os{opath, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(seqan3::interleaved_bloom_filter(ibf));
}
