#pragma once

#include <deque>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/view/all.hpp>
#include <seqan3/range/view/kmer_hash.hpp>
#include <seqan3/search/kmer_index/shape.hpp>
#include <seqan3/std/ranges>

// k-mer size
uint8_t p_k{20};
// window size
uint16_t p_w{60};
// shape
uint64_t p_shape = 0;
// Random, but static value for xor for hashes. Counteracts consecutive minimizers.
// E.g., without it, the next minimizer after a poly-A region AAAAA would be most likely something like AAAAC.
// Setting seed to 0, will lead to lexicographically smallest k-mer
uint64_t p_seed{0x8F3F73B5CF1C9ADE};

struct cmd_arguments
{
    // ibf specific arguments
    std::vector<std::filesystem::path> sequence_files;
    std::filesystem::path genome_file;
    std::vector<size_t> bits{};
    size_t num_hash{1};
    std::filesystem::path path_out{"./"};
    std::vector<float> expression_levels{}; // 0.5,1,2,4
    std::vector<int> samples{};
    std::string aggregate_by{"median"};
    size_t random{10};

    // search specific arguments
    std::filesystem::path gene_file;
    std::filesystem::path exp_file;
    std::filesystem::path path_in{"./"};
    float expression{1.0};

    // Arguments used for ibf and search
    bool compressed = false;
    uint8_t k{20};
    uint16_t window_size{60};
    uint64_t shape;
    uint64_t seed{0x8F3F73B5CF1C9ADE};

};

struct my_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;               // instead of dna5
};

struct Minimizer
{
public:
    // k-mer size
    uint8_t k = p_k;
    // window size
    uint16_t w = p_w;
    // shape
    seqan3::shape shape{seqan3::ungapped{k}};
    // seed
    uint64_t seed = p_seed;
    // start positions of minimizers
    std::vector<uint64_t> minBegin;
    // end positions of minimizers
    std::vector<uint64_t> minEnd;


    inline void resize(uint16_t newkmerPos, uint32_t neww, uint64_t newShape = p_shape, uint64_t newSeed = p_seed)
    {
        k = newkmerPos;
        w = neww;
        if (newShape == 0)
            shape = seqan3::ungapped{k};
        else
            shape = seqan3::bin_literal{newShape};
        seed = newSeed;
    }

    std::vector<uint64_t> getMinimizer(seqan3::dna4_vector const & text)
    {
        if (k > text.size())
            return std::vector<uint64_t> {};

        // Reverse complement without copying/modifying the original string
        seqan3::dna4_vector revComp = text | std::view::reverse | seqan3::view::complement;

        uint64_t possible = text.size() > w ? text.size() - w + 1 : 1;
        uint32_t windowKmers = w - k + 1;
        uint64_t kmerPos = text.size() - k + 1;

        std::vector<uint64_t> kmerHashes{};
        kmerHashes.reserve(possible); // maybe rather reserve to expected?

        // Stores hash, begin and end for all k-mers in the window
        std::deque<uint64_t> windowValues;

        auto kmerHashing = text | seqan3::view::kmer_hash(shape);
        auto revcHashing = revComp | seqan3::view::kmer_hash(shape) | std::view::reverse;

        uint32_t i = 0;

        // Initialisation. We need to compute all hashes for the first window.
        for ( auto && [hf, hr] : std::ranges::view::zip(kmerHashing | seqan3::view::take_exactly(windowKmers),
             revcHashing | seqan3::view::take_exactly(windowKmers)) )
        {
            // Get smallest canonical k-mer
            uint64_t kmerHash = hf ^ seed;
            uint64_t revcHash = hr ^ seed;
            if (kmerHash <= revcHash)
            {
                windowValues.push_back(kmerHash);

            }
            else
            {
                windowValues.push_back(revcHash);
            }
        }

        auto min = std::min_element(std::begin(windowValues), std::end(windowValues));
        kmerHashes.push_back(*min);

        // For the following windows, we remove the first window k-mer (is now not in window) and add the new k-mer
        // that results from the window shifting
        bool minimizer_changed{false};
        for ( auto && [hf, hr] : std::ranges::view::zip(kmerHashing | seqan3::view::drop(windowKmers),
             revcHashing | seqan3::view::drop(windowKmers)))
        {
            if (min == std::begin(windowValues))
            {
                windowValues.pop_front();
                min = std::min_element(std::begin(windowValues), std::end(windowValues));
                minimizer_changed = true;
            }
            else
                windowValues.pop_front();

            uint64_t kmerHash = hf ^ seed;
            uint64_t revcHash = hr ^ seed;
            if (kmerHash <= revcHash)
                windowValues.push_back(kmerHash);
            else
                windowValues.push_back(revcHash);

            if (windowValues.back() < *min)
            {
                min = std::end(windowValues) - 1;
                minimizer_changed = true;
            }

            if (minimizer_changed)
            {
                kmerHashes.push_back(*min);
                minimizer_changed = false;
            }
        }

        return kmerHashes;
    }
};

std::vector<uint64_t> compute_minimizer(const seqan3::dna4_vector & seq, uint8_t k = p_k, uint16_t window_size = p_w,
                                        uint64_t shape = p_shape, uint64_t seed = p_seed)
{
//  size_t const window_size = ((length(seq) - k) / LeastNumMinimizers) + k; // ensure at least 3 minimizer
    if (k > window_size)
        throw std::logic_error("Cannot divide the sequence into the specified window size."); // return std::array<uint64_t, 3> {};

    Minimizer minimizer;
    minimizer.resize(k, window_size, shape, seed);

    return minimizer.getMinimizer(seq);
}

auto compute_occurrences(const std::vector<seqan3::dna4_vector> & seqs, uint8_t k = p_k, uint16_t window_size = p_w,
                                        uint64_t shape = p_shape, uint64_t seed = p_seed)
{
    std::unordered_map<uint64_t, uint32_t> occurring_kmers{};
    unsigned total_count{0};

    for (auto & seq : seqs)
    {
        for (auto & minHash : compute_minimizer(seq, k, window_size,shape, seed))
        {
            ++occurring_kmers[minHash];
        }
    }

    return occurring_kmers;
}
