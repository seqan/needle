// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/needle/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <robin_hood.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/search/kmer_index/shape.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

inline constexpr static uint64_t adjust_seed(uint8_t const kmer_size, uint64_t const seed = 0x8F3F73B5CF1C9ADEULL) noexcept
{
    return seed >> (64u - 2u * kmer_size);
}

//!\brief arguments used for all tools
struct all_arguments
{
    std::filesystem::path path_out{"./"};
    uint8_t threads{1};
};

//!\brief arguments used for estimate, ibf, minimiser
struct min_arguments : all_arguments
{
    uint8_t k{20};
    seqan3::seed s{0x8F3F73B5CF1C9ADEULL};
    seqan3::shape shape = seqan3::ungapped{k};
    seqan3::window_size w_size{60};
};

//!\brief arguments used for estimate, ibf, ibfmin
struct estimate_ibf_arguments : min_arguments
{
    bool compressed = false;
    std::vector<uint16_t> expression_thresholds{}; // Expression levels which should be created
    uint8_t number_expression_thresholds{}; // If set, the expression levels are determined by the program.
    bool samplewise{false};

    template<class Archive>
    void save(Archive & archive) const
    {
        archive(k);
        archive(w_size.get());
        archive(s.get());
        archive(shape);
        archive(compressed);
        archive(number_expression_thresholds);
        archive(expression_thresholds);
        archive(samplewise);
    }

    template<class Archive>
    void load(Archive & archive)
    {
        archive(k);
        archive(w_size.get());
        archive(s.get());
        archive(shape);
        archive(compressed);
        archive(number_expression_thresholds);
        archive(expression_thresholds);
        archive(samplewise);
    }
};

/*! \brief Function, loading arguments
 *  \param args   arguments to load
 *  \param ipath Path, where the arguments can be found.
 */
[[maybe_unused]] static void load_args(estimate_ibf_arguments & args, std::filesystem::path ipath)
{
    std::ifstream is{ipath, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    iarchive(args);
}

/*! \brief Function, which stores the arguments
 *  \param args  arguments to store
 *  \param opath Path, where the arguments should be stored.
 */
[[maybe_unused]] static void store_args(estimate_ibf_arguments const & args, std::filesystem::path opath)
{
    std::ofstream os{opath, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(args);
}

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
 *  \param opath Path, where the IBF should be stored.
 */
template <class IBFType>
void store_ibf(IBFType const & ibf,
               std::filesystem::path opath)
{
    std::ofstream os{opath, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(ibf);
}
