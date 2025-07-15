// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/kmer_index/shape.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

#include <hibf/hierarchical_interleaved_bloom_filter.hpp>

static inline constexpr uint64_t adjust_seed(uint8_t const kmer_size,
                                             uint64_t const seed = 0x8F'3F'73'B5'CF'1C'9A'DEULL) noexcept
{
    return seed >> (64u - 2u * kmer_size);
}

//!\brief arguments used for all tools
struct all_arguments
{
    std::filesystem::path path_out{"./"};
    uint16_t threads{1u};
};

//!\brief arguments used for estimate, ibf, minimiser
struct minimiser_arguments : all_arguments
{
    uint8_t k{20};
    seqan3::seed s{0x8F'3F'73'B5'CF'1C'9A'DEULL};
    seqan3::shape shape = seqan3::ungapped{k};
    seqan3::window_size w_size{60};
};

//!\brief arguments used for estimate, ibf, ibfmin
struct estimate_ibf_arguments : minimiser_arguments
{
    bool compressed = false;
    std::vector<uint16_t> expression_thresholds{}; // Expression levels which should be created
    uint8_t number_expression_thresholds{};        // If set, the expression levels are determined by the program.
    bool samplewise{false};

    template <typename archive_t>
    void serialize(archive_t & archive)
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

struct minimiser_file_input_arguments
{
    // Needs to be defined when only minimisers appearing in this file should be stored
    std::filesystem::path include_file;
    // Needs to be defined when minimisers appearing in this file should NOT be stored
    std::filesystem::path exclude_file;
    std::vector<size_t> samples{}; // Can be used to indicate that sequence files belong to the same experiment
    bool paired = false;           // If true, than experiments are seen as paired-end experiments
    bool experiment_names = false; // Flag, if names of experiment should be stored in a txt file
    bool ram_friendly = false;
};

/*! \brief Function, loading arguments
 *  \param args   arguments to load
 *  \param ipath Path, where the arguments can be found.
 */
inline void load_args(estimate_ibf_arguments & args, std::filesystem::path const & ipath)
{
    std::ifstream is{ipath, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    iarchive(args);
}

/*! \brief Function, which stores the arguments
 *  \param args  arguments to store
 *  \param opath Path, where the arguments should be stored.
 */
inline void store_args(estimate_ibf_arguments const & args, std::filesystem::path const & opath)
{
    std::ofstream os{opath, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(args);
}

//!\brief Use dna4 instead of default dna5
struct dna4_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
    //TODO: Should I use a bitcompressed_vector to save memory but with the disadvantage of losing speed?
    //template <typename alph>
    //using sequence_container = seqan3::bitcompressed_vector<alph>;
};

using sequence_file_t = seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>>;
using sequence_file_with_id_t =
    seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::id, seqan3::field::seq>>;

/*! \brief Function, loading compressed and uncompressed ibfs
 *  \param ibf   ibf to load
 *  \param ipath Path, where the ibf can be found.
 */
template <typename ibf_t>
void load_ibf(ibf_t & ibf, std::filesystem::path const & ipath)
{
    std::ifstream is{ipath, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    iarchive(ibf);
}

/*! \brief Function, which stored compressed and uncompressed ibfs
 *  \param ibf   The IBF to store.
 *  \param opath Path, where the IBF should be stored.
 */
template <typename ibf_t>
void store_ibf(ibf_t const & ibf, std::filesystem::path const & opath)
{
    std::ofstream os{opath, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(ibf);
}
