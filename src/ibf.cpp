// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/needle/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <chrono>
#include <deque>
#include <iostream>
#include <math.h>
#include <mutex>
#include <numeric>
#include <omp.h>
#include <string>
#include <thread>
#include <algorithm> //reorded because of this error:https://github.com/Homebrew/homebrew-core/issues/44579

#include <filesystem>
#include <ranges>

#if SEQAN3_WITH_CEREAL
#include <cereal/archives/binary.hpp>
#endif // SEQAN3_WITH_CEREAL

#include <robin_hood.h>

#include <seqan3/alphabet/container/concatenated_sequences.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/contrib/parallel/buffer_queue.hpp>
#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/io/stream/detail/fast_istreambuf_iterator.hpp>
#include <seqan3/utility/container/dynamic_bitset.hpp>
#include <seqan3/utility/parallel/detail/latch.hpp>

#include "ibf.h"
#include "shared.h"

// Create set with hashes from the minimisers from an include or exclude file.
void get_include_set_table(min_arguments const & args, std::filesystem::path const include_file,
                           robin_hood::unordered_set<uint64_t> & include_table)
{
    seqan3::sequence_file_input<my_traits,  seqan3::fields<seqan3::field::seq>> fin3{include_file};
    for (auto & [seq] : fin3)
    {
        if (seq.size() >= args.w_size.get())
        {
            for (auto && minHash : seqan3::views::minimiser_hash(seq, args.shape, args.w_size, args.s))
                include_table.insert(minHash);
        }
    }
}

// Chech if file has fasta format to estimate cutoffs.
inline bool check_for_fasta_format(std::vector<std::string> const & valid_extensions, std::string const & file_path)
{

    auto case_insensitive_string_ends_with = [&] (std::string_view str, std::string_view suffix)
    {
        size_t const suffix_length{suffix.size()};
        size_t const str_length{str.size()};
        return suffix_length > str_length ?
               false :
               std::ranges::equal(str.substr(str_length - suffix_length), suffix, [] (char const chr1, char const chr2)
               {
                   return std::tolower(chr1) == std::tolower(chr2);
               });
    };

    auto case_insensitive_ends_with = [&] (std::string const & ext)
    {
        return case_insensitive_string_ends_with(file_path, ext);
    };

    return std::ranges::find_if(valid_extensions, case_insensitive_ends_with) != valid_extensions.end();
}

// Determine cutoff for one experiment
uint8_t calculate_cutoff(std::filesystem::path sequence_file, int samples)
{
    // Cutoff according to Mantis paper -1 because we use "<" and not "<="
    uint16_t const default_cutoff{49};
    uint8_t cutoff{default_cutoff};
    std::array<uint16_t, 4> const cutoffs{0, 2, 9, 19};
    std::array<uint64_t, 4> const cutoff_bounds{314'572'800, 524'288'000, 1'073'741'824, 3'221'225'472};
    cutoff = default_cutoff;

    // Since the curoffs are based on the filesize of a gzipped fastq file, we try account for the other cases:
    // We multiply by two if we have fasta input.
    // We divide by 3 if the input is not compressed.
    bool const is_compressed = sequence_file.extension() == ".gz" || sequence_file.extension() == ".bgzf" || sequence_file.extension() == ".bz2";
    bool const is_fasta = is_compressed ? check_for_fasta_format(seqan3::format_fasta::file_extensions, sequence_file.stem())
                                       : check_for_fasta_format(seqan3::format_fasta::file_extensions, sequence_file.extension());
    size_t const filesize = std::filesystem::file_size(sequence_file) * samples * (is_fasta ? 2 : 1) / (is_compressed ? 1 : 3);

    for (size_t k = 0; k < cutoff_bounds.size(); ++k)
    {
        if (filesize <= cutoff_bounds[k])
        {
            cutoff = cutoffs[k];
            break;
        }
    }
    return cutoff;
}

// Fill hash table with minimisers greater than the cutoff.
void fill_hash_table(min_arguments const & args,
                     seqan3::sequence_file_input<my_traits,  seqan3::fields<seqan3::field::seq>> & fin,
                     robin_hood::unordered_node_map<uint64_t, uint16_t> & hash_table,
                     robin_hood::unordered_node_map<uint64_t, uint8_t> & cutoff_table,
                     robin_hood::unordered_set<uint64_t> const & include_set_table,
                     robin_hood::unordered_set<uint64_t> const & exclude_set_table,
                     bool const only_include = false, uint8_t cutoff = 0)
{
    for (auto & [seq] : fin)
    {
        for (auto && minHash : seqan3::views::minimiser_hash(seq, args.shape, args.w_size, args.s))
        {
            if ((only_include & (include_set_table.contains(minHash))) | (!only_include) & !(exclude_set_table.contains(minHash)))
            {
                auto it = hash_table.find(minHash);
                // If minHash is already in hash table, increase count in hash table
                if (it != hash_table.end())
                {
                    it->second = std::min<uint16_t>(65534u, hash_table[minHash] + 1);
                }
                // If minHash equals now the cutoff than add it to the hash table and add plus one for the current
                // iteration.
                else if (cutoff_table[minHash] == cutoff)
                {
                    hash_table[minHash] = cutoff_table[minHash] + 1;
                    cutoff_table.erase(minHash);
                }
                else if (cutoff == 0)
                {
                    hash_table[minHash]++;
                }
                // If none of the above, increase count in cutoff table. Cutoff Table increases RAM usage by storing
                // minimisers with a low occurence in a smaller hash table.
                else
                {
                    cutoff_table[minHash]++;
                }
            }
        }
    }
}

void fill_hash_table_parallel(min_arguments const & args,
                              seqan3::sequence_file_input<my_traits,  seqan3::fields<seqan3::field::seq>> & fin,
                              robin_hood::unordered_node_map<uint64_t, uint16_t> & hash_table,
                              robin_hood::unordered_node_map<uint64_t, uint8_t> & cutoff_table,
                              robin_hood::unordered_set<uint64_t> const & include_set_table,
                              robin_hood::unordered_set<uint64_t> const & exclude_set_table,
                              bool const only_include = false,
                              uint8_t cutoff = 0)
{
    using value_t = typename robin_hood::unordered_node_map<uint64_t, uint16_t>::value_type;
    using value_t2 = typename robin_hood::unordered_node_map<uint64_t, uint8_t>::value_type;
    using local_hash_table_t = robin_hood::unordered_node_map<uint64_t, std::atomic<uint16_t>>;
    using local_cutoff_table_t = robin_hood::unordered_node_map<uint64_t, std::atomic<uint8_t>>;

    // Step 1: load file in batches
    // Step 2: construct minimiser hash values
    // Step 3: Sort minimiser hash values but pertain their original positions
    // Step 4: Provid simple table with minimiser value and count pair
    // Step 5: Store all non-referenced minimiser values per batch

    size_t thread_count = args.threads;
    // Run thread: block execution and load next 10000 sequences.

    auto seq_file_it = std::ranges::begin(fin);
    size_t chunk_count{};
    using sequence_t = seqan3::dna4_vector;

    auto load_next_chunk = [&] ()
    {
        constexpr size_t batch_size = 100000;

        std::vector<sequence_t> sequence_batch{};
        sequence_batch.reserve(batch_size);

        while (seq_file_it != std::ranges::end(fin) && sequence_batch.size() < batch_size)
        {
            sequence_batch.push_back((*seq_file_it).sequence());
            ++seq_file_it;
        }

        return sequence_batch;
    };

    auto count_minimiser = [&] (auto & local_hash_table, std::vector<uint64_t> minimisers)
    {
        // Sort the minimiser by their value.
        std::ranges::sort(minimisers, std::less{});

        // Fill the table with all minimiser counts.
        auto minimiser_it = minimisers.begin();
        auto minimiser_end = minimisers.end();

        std::vector<uint64_t> orphaned_minimiser{};

        while (minimiser_it != minimiser_end)
        {
            uint64_t current_minimiser = *minimiser_it;
            auto predicate = [=] (uint64_t const other_hash) { return other_hash == current_minimiser; };
            auto next_minimiser_it = std::ranges::find_if_not(minimiser_it, minimiser_end, predicate);
            // minimiser_it now points to the first non equal position
            size_t const minimiser_count = std::ranges::distance(minimiser_it, next_minimiser_it);

            if ((only_include && (include_set_table.contains(current_minimiser))) || (!only_include) && !(exclude_set_table.contains(current_minimiser)))
            {
                if (auto it = local_hash_table.find(current_minimiser); it != local_hash_table.end()) // update
                {
                    it->second = static_cast<uint16_t>(std::min(65534ul, it->second + minimiser_count));
                }
                else if (minimiser_count > cutoff)
                {
                     // insert first.
                    local_hash_table[current_minimiser] = minimiser_count;
                }
                else // not above cutoff.
                {
                    orphaned_minimiser.insert(orphaned_minimiser.end(), minimiser_it, next_minimiser_it);
                }
            }
            minimiser_it = next_minimiser_it;
        }
        return orphaned_minimiser;
    };

    // Block 1:
    std::vector<std::vector<uint64_t>> thread_local_remaining_minimisers{};
    thread_local_remaining_minimisers.resize(thread_count);
    std::vector<local_hash_table_t> thread_local_hash_tables{};
    thread_local_hash_tables.resize(thread_count);
    bool is_merged{false};
    bool fill_queue{true};

    seqan3::detail::latch sync_point{thread_count};
    seqan3::detail::latch sync_point_2{thread_count};
    std::atomic<size_t> remaining_minimisers_size{};
    std::mutex load_mutex{};
    std::vector<std::pair<size_t, size_t>> intervals{};
    std::optional<seqan3::contrib::fixed_buffer_queue<std::pair<size_t, size_t>>> queue;

    auto job = [&] (size_t const thread_id)
    {
        while (true)
        {
            std::vector<sequence_t> sequence_batch{};
            { // critical region
                std::scoped_lock load_lk{load_mutex};
                sequence_batch = load_next_chunk();
            }
            if (sequence_batch.empty()) // Stop construction: no more elements are coming
                break;

            // Construct the set of all minimisers for all sequences.
            std::vector<uint64_t> minimisers{};
            minimisers.reserve(sequence_batch.size() * (sequence_batch[0].size() - args.w_size.get() + 1));
            for (auto & sequence : sequence_batch)
                std::ranges::move(sequence | seqan3::views::minimiser_hash(args.shape, args.w_size, args.s),
                                  std::back_inserter(minimisers));

            auto orphaned_minimiser = count_minimiser(thread_local_hash_tables[thread_id], std::move(minimisers));
            thread_local_remaining_minimisers[thread_id].insert(thread_local_remaining_minimisers[thread_id].end(),
                                                                orphaned_minimiser.begin(),
                                                                orphaned_minimiser.end());
        }

        sync_point.arrive_and_wait();

        {// sequential phase to merge sub tables.
            std::scoped_lock lk{load_mutex};
            if (!is_merged)
            {
                // Merge local hash_tables.
                for (auto & local_hash_table : thread_local_hash_tables)
                {
                    for (auto && [key, counter] : local_hash_table)
                    {
                        if (auto it = hash_table.find(key); it != hash_table.end())
                            it->second = static_cast<uint16_t>(std::min<uint32_t>(65534ul, it->second + counter.load()));
                        else
                            hash_table.insert(value_t{key, counter.load()});
                    }

                    local_hash_table.clear();
                }

                std::vector<uint64_t> remaining_minimisers{};
                for (auto & local_remaining_minimisers : thread_local_remaining_minimisers)
                {
                    std::vector<uint64_t> local_remaining_minimisers2 =
                    count_minimiser(hash_table, std::move(local_remaining_minimisers));
                    remaining_minimisers.insert(remaining_minimisers.end(), local_remaining_minimisers2.begin(), local_remaining_minimisers2.end());
                }

                std::ranges::sort(remaining_minimisers, std::less{});
                auto minimiser_it = remaining_minimisers.begin();
                auto minimiser_end = remaining_minimisers.end();

                while (minimiser_it != minimiser_end)
                {
                    uint64_t current_minimiser = *minimiser_it;
                    auto predicate = [=] (uint64_t const other_hash) { return other_hash == current_minimiser; };
                    auto next_minimiser_it = std::ranges::find_if_not(minimiser_it, minimiser_end, predicate);
                    // minimiser_it now points to the first non equal position
                    size_t const minimiser_count = std::ranges::distance(minimiser_it, next_minimiser_it);

                    if ((only_include && (include_set_table.contains(current_minimiser))) || (!only_include) && !(exclude_set_table.contains(current_minimiser)))
                    {
                        if (auto it = hash_table.find(current_minimiser); it != hash_table.end()) // update
                        {
                            it->second = static_cast<uint16_t>(std::min(65534ul, it->second + minimiser_count));
                        }
                        else if (minimiser_count > cutoff)
                        {
                             // insert first.
                            hash_table[current_minimiser] = minimiser_count;
                        }
                        else if (auto it = cutoff_table.find(current_minimiser); it != cutoff_table.end())
                        {
                            if ((it->second + minimiser_count) > cutoff)
                            {
                                hash_table[current_minimiser] = cutoff_table[current_minimiser] + minimiser_count;
                                cutoff_table.erase(current_minimiser);
                            }
                            else
                            {
                                it->second = it->second + minimiser_count;
                            }
                        }
                        else
                        {
                            cutoff_table[current_minimiser] = minimiser_count;
                        }
                    }
                    minimiser_it = next_minimiser_it;
                }
                is_merged = true;
            }
        }
    };

    std::vector<std::thread> thread_pool{};
    for (size_t i = 0; i < thread_count; ++i)
        thread_pool.emplace_back(job, i);

    // Wait for all threads to finish.
    for (auto & thread : thread_pool)
        if (thread.joinable())
            thread.join();
}

void count_genome(min_arguments const & args, std::filesystem::path include_file,
                  std::filesystem::path exclude_file)
{
    robin_hood::unordered_set<uint64_t> include_set_table{};
    robin_hood::unordered_set<uint64_t> exclude_set_table{};
    std::ofstream outfile;

    if (exclude_file != "")
    {
        seqan3::sequence_file_input<my_traits,  seqan3::fields<seqan3::field::seq>> fin{exclude_file};
        for (auto & [seq] : fin)
        {
            if (seq.size() >= args.w_size.get())
            {
                for (auto && minHash : seqan3::views::minimiser_hash(seq, args.shape, args.w_size, args.s))
                    exclude_set_table.insert(minHash);
            }
        }
    }

    seqan3::sequence_file_input<my_traits,  seqan3::fields<seqan3::field::seq>> fin2{include_file};
    for (auto & [seq] : fin2)
    {
        if (seq.size() >= args.w_size.get())
        {
            for (auto && minHash : seqan3::views::minimiser_hash(seq, args.shape, args.w_size, args.s))
            {
                if ( !(exclude_set_table.contains(minHash)))
                    include_set_table.insert(minHash);
            }
        }
    }

    // Write minimiser to binary
    outfile.open(std::string{args.path_out} + std::string{include_file.stem()}
                 + ".genome", std::ios::binary);

    for (auto && hash : include_set_table)
    {
        outfile.write(reinterpret_cast<const char*>(&hash), sizeof(hash));
    }
    outfile.close();
}

void count(min_arguments const & args, std::vector<std::filesystem::path> sequence_files, std::filesystem::path include_file,
           std::filesystem::path genome_file, bool paired)
{
    robin_hood::unordered_node_map<uint64_t, uint16_t> hash_table{};
    // Create a smaller cutoff table to save RAM, this cutoff table is only used for constructing the hash table
    // and afterwards discarded.
    robin_hood::unordered_node_map<uint64_t, uint8_t>  cutoff_table;
    robin_hood::unordered_set<uint64_t> include_set_table{};
    robin_hood::unordered_set<uint64_t> exclude_set_table{};
    std::vector<uint64_t> counter{};
    uint64_t exp{};
    std::ifstream infile;
    std::ofstream outfile;
    int j;

    // Read minimiser from binary
    infile.open(genome_file, std::ios::binary);
    uint64_t minimiser;
    while(infile.read((char*)&minimiser, sizeof(minimiser)))
    {
        include_set_table.insert(minimiser);
    }
    infile.close();

    for (unsigned i = 0; i < sequence_files.size(); i++)
    {

        if (paired)
        {
            seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{sequence_files[i]};
            fill_hash_table(args, fin, hash_table, cutoff_table, include_set_table, exclude_set_table, true);
            i++;
            fin = sequence_files[i];
            fill_hash_table(args, fin, hash_table, cutoff_table, include_set_table, exclude_set_table, true);
        }
        else
        {
            seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{sequence_files[i]};
            fill_hash_table(args, fin, hash_table, cutoff_table, include_set_table, exclude_set_table, true);
        }
        cutoff_table.clear();

        outfile.open(std::string{args.path_out} + std::string{sequence_files[i].stem()} + ".count.out");
        j = 0;
        seqan3::sequence_file_input<my_traits,  seqan3::fields<seqan3::field::id, seqan3::field::seq>> fin2{include_file};
        for (auto & [id, seq] : fin2)
        {
            if (seq.size() >= args.w_size.get())
            {
                for (auto && minHash : seqan3::views::minimiser_hash(seq, args.shape, args.w_size, args.s))
                    counter.push_back(hash_table[minHash]);
                std::nth_element(counter.begin(), counter.begin() + counter.size()/2, counter.end());
                exp =  counter[counter.size()/2];
                counter.clear();
                outfile << id << "\t" << exp << "\n";
                ++j;
            }
        }
        outfile.close();
        hash_table.clear();
    }
}

void read_binary(std::filesystem::path filename, robin_hood::unordered_node_map<uint64_t, uint16_t> & hash_table)
{
    std::ifstream fin;

    uint8_t small_buffer;
    uint32_t window;
    uint64_t buffer;
    fin.open(filename, std::ios::binary);
    fin.read((char*)&buffer, sizeof(buffer));
    fin.read((char*)&small_buffer, sizeof(small_buffer));
    fin.read((char*)&small_buffer, sizeof(small_buffer));
    fin.read((char*)&window, sizeof(window));
    fin.read((char*)&buffer, sizeof(buffer));
    bool ungapped;
    fin.read((char*)&ungapped, sizeof(ungapped));
    if (!ungapped)
    {
        fin.read((char*)&buffer, sizeof(buffer));
    }

    uint64_t minimiser;
    uint16_t minimiser_count;

    while(fin.read((char*)&minimiser, sizeof(minimiser)))
    {
        fin.read((char*)&minimiser_count, sizeof(minimiser_count));
        hash_table[minimiser] = minimiser_count;
    }

    fin.close();
}

void read_binary_start(min_arguments & args,
                 std::filesystem::path filename,
                 uint64_t & num_of_minimisers, uint8_t & cutoff)
{
    std::ifstream fin;

    uint32_t window;
    uint64_t buffer;
    uint8_t small_buffer;
    fin.open(filename, std::ios::binary);
    fin.read((char*)&buffer, sizeof(buffer));
    num_of_minimisers = buffer;

    fin.read((char*)&small_buffer, sizeof(small_buffer));
    cutoff = small_buffer;
    fin.read((char*)&args.k, sizeof(args.k));
    fin.read((char*)&window, sizeof(window));
    args.w_size = seqan3::window_size{window};
    fin.read((char*)&buffer, sizeof(buffer));
    args.s = seqan3::seed{buffer};

    bool ungapped;
    fin.read((char*)&ungapped, sizeof(ungapped));
    if (ungapped)
    {
        args.shape = seqan3::ungapped{args.k};
    }
    else
    {
        fin.read((char*)&buffer, sizeof(buffer));
        args.shape = seqan3::bin_literal{buffer};
    }

    fin.close();
}

// Check number of expression levels, sort expression levels
void check_expression(std::vector<uint16_t> & expression_thresholds, uint8_t & number_expression_thresholds,
                      std::filesystem::path const expression_by_genome_file)
{
    // Sort given expression rates
    sort(expression_thresholds.begin(), expression_thresholds.end());

     // If no expression levels are given and the no number of expression levels is specified, throw.
    if ((number_expression_thresholds == 0) & (expression_thresholds.size() == 0))
    {
        throw std::invalid_argument{"Error. Please set the expression levels OR give the number of expression levels."};
    }
    else if ((expression_by_genome_file != "") & (expression_thresholds.size() > 0))
    {
        throw std::invalid_argument{"Error. The determination of expression levels can not be used with individual levels"
                                    " already given. Please set the expression levels without the option "
                                    "--level-by-genome OR use the number of expression levels with that option."};
    }
    else if (number_expression_thresholds == 0)
    {
        number_expression_thresholds = expression_thresholds.size();
    }
    else if ((number_expression_thresholds != expression_thresholds.size()) & (expression_thresholds.size() > 0))
    {
        throw std::invalid_argument{"Error. Please set the expression levels OR give the number of expression levels."};
    }

}

// Check and set samples and cutoffs
void check_cutoffs_samples(std::vector<std::filesystem::path> const & sequence_files,
                           bool const paired, std::vector<int> & samples,
                           std::vector<uint8_t> & cutoffs)
{
    if (paired) // If paired is true, a pair is seen as one sample
        samples.assign(sequence_files.size()/2,2);
    if (samples.empty()) // If no samples are given and not paired, every file is seen as one experiment
        samples.assign(sequence_files.size(),1);
    if (cutoffs.size() == 1) // If one cutoff is given, every experiment gets this cutoff.
        cutoffs.assign(samples.size(), cutoffs[0]);

    // If sum of minimiser_args.samples is not equal to number of files, throw error
    else if (std::accumulate(samples.rbegin(), samples.rend(), 0) != sequence_files.size())
        throw std::invalid_argument{"Error. Incorrect command line input for multiple-samples."};
}

// Check input of fpr
void check_fpr(uint8_t const number_expression_thresholds, std::vector<double> & fprs)
{
    // If no bin size is given or not the right amount, throw error.
    if (fprs.empty())
    {
        throw std::invalid_argument{"Error. Please give a false positive rate for the IBFs."};
    }
    // If only one ibf size is given, set it for all thresholds.
    if (fprs.size() == 1)
    {
        fprs.assign(number_expression_thresholds, fprs[0]);
    }
    else if (fprs.size() != number_expression_thresholds)
    {
        throw std::invalid_argument{"Error. Length of false positive rates for IBFs is not equal to length of expression "
                                    "thresholds."};
    }
}

// Calculate expression thresholds and sizes
void get_expression_thresholds(uint8_t const number_expression_thresholds,
                           robin_hood::unordered_node_map<uint64_t, uint16_t> const & hash_table,
                           std::vector<uint16_t> & expression_thresholds, std::vector<uint64_t> & sizes,
                           robin_hood::unordered_set<uint64_t> const & genome, uint8_t cutoff, bool all = true)
{
    // Calculate expression thresholds by taking median recursively
    std::vector<uint16_t> counts;
    for (auto && elem : hash_table)
    {
        if (all | genome.contains(elem.first))
            counts.push_back(elem.second);
    }

    std::size_t dev{2};
    std::size_t prev_pos{0};
    auto prev_exp{0};
    auto exp{0};
    auto max_elem = *std::max_element(counts.begin(), counts.end());
    // Zero Level = cutoff + 1
    expression_thresholds.push_back(cutoff + 1);
    // First Level
    std::nth_element(counts.begin() + prev_pos, counts.begin() +  prev_pos + counts.size()/dev, counts.end());
    exp = counts[prev_pos + counts.size()/dev];
    prev_pos = prev_pos + counts.size()/dev;
    dev = dev*2;
    expression_thresholds.push_back(exp);
    sizes.push_back(prev_pos);

    while((expression_thresholds.size() < number_expression_thresholds) & (prev_exp < max_elem) & (dev < counts.size()))
    {
        std::nth_element(counts.begin() + prev_pos, counts.begin() +  prev_pos + counts.size()/dev, counts.end());
        exp = counts[prev_pos + counts.size()/dev];
        prev_pos = prev_pos + counts.size()/dev;
        dev = dev*2;

        // If expression does not change compared to previous one, do not store it again as an expression threshold.
        if ((exp - prev_exp) > 1)
        {
            expression_thresholds.push_back(exp);
            sizes.push_back(prev_pos);
        }

        prev_exp = exp;
    }
    sizes.push_back(prev_pos);
    // In case not all levels have a threshold, give the last levels a maximal threshold, which can not be met by any minimiser.
    while(expression_thresholds.size() < number_expression_thresholds)
        expression_thresholds.push_back(max_elem + 1);
    counts.clear();
}

// Estimate the file size for every expression level, necessary when samplewise=false, because then it is completly
// unclear how many minimisers are to store per file.
void get_filsize_per_expression_level(std::filesystem::path filename, uint8_t const number_expression_thresholds,
                                      std::vector<uint16_t> const & expression_thresholds, std::vector<uint64_t> & sizes,
                                      robin_hood::unordered_set<uint64_t> const & genome, bool all = true)
{
    std::ifstream fin;
    uint8_t small_buffer;
    uint32_t window;
    uint64_t buffer;
    fin.open(filename, std::ios::binary);
    fin.read((char*)&buffer, sizeof(buffer));
    fin.read((char*)&small_buffer, sizeof(small_buffer));
    fin.read((char*)&small_buffer, sizeof(small_buffer));
    fin.read((char*)&window, sizeof(window));
    fin.read((char*)&buffer, sizeof(buffer));
    bool ungapped;
    fin.read((char*)&ungapped, sizeof(ungapped));
    if (!ungapped)
    {
        fin.read((char*)&buffer, sizeof(buffer));
    }

    uint64_t minimiser;
    uint16_t minimiser_count;
    sizes.assign(number_expression_thresholds, 0);

    while(fin.read((char*)&minimiser, sizeof(minimiser)))
    {
        fin.read((char*)&minimiser_count, sizeof(minimiser_count));
        if (all | genome.contains(minimiser))
        {
            // Find the level with the smallest greater value than the minimiser occurrence, in the level before that the
            // minimiser is going to be stored.
            auto p = std::upper_bound(expression_thresholds.begin(), expression_thresholds.end(), minimiser_count);
            if(p != expression_thresholds.begin())
                sizes[(p-expression_thresholds.begin())-1]++;
        }
    }

    fin.close();
}

// Actual ibf construction
template<bool samplewise, bool minimiser_files_given = true>
void ibf_helper(std::vector<std::filesystem::path> const & minimiser_files,
                std::vector<double> const & fprs,
                estimate_ibf_arguments & ibf_args, std::vector<uint8_t> & cutoffs = {},
                size_t num_hash = 1, std::filesystem::path expression_by_genome_file = "",
                minimiser_arguments const & minimiser_args = {})
{

    size_t num_files;
    if constexpr (minimiser_files_given)
        num_files = minimiser_files.size();
    else
        num_files = minimiser_args.samples.size();

    std::vector<std::vector<uint16_t>> expressions{};
    std::vector<std::vector<uint64_t>> sizes{};
    sizes.assign(num_files, {});

    bool const calculate_cutoffs = cutoffs.empty();

    robin_hood::unordered_set<uint64_t> include_set_table; // Storage for minimisers in include file
    robin_hood::unordered_set<uint64_t> exclude_set_table; // Storage for minimisers in exclude file
    if constexpr(samplewise)
    {
        std::vector<uint16_t> zero_vector(ibf_args.number_expression_thresholds);
        for (unsigned j = 0; j < num_files; j++)
            expressions.push_back(zero_vector);

    }
    if constexpr (!minimiser_files_given)
    {
        if (minimiser_args.include_file != "")
            get_include_set_table(ibf_args, minimiser_args.include_file, include_set_table);
        if (minimiser_args.exclude_file != "")
            get_include_set_table(ibf_args, minimiser_args.exclude_file, exclude_set_table);
    }

    if (minimiser_args.ram_friendly)
        omp_set_num_threads(1);
    else
        omp_set_num_threads(ibf_args.threads);
    seqan3::contrib::bgzf_thread_count = ibf_args.threads;

    size_t const chunk_size = std::clamp<size_t>(std::bit_ceil(num_files / ibf_args.threads), 8u, 64u);

    // If expression_thresholds should only be depending on minimsers in a certain genome file, genome is created.
    robin_hood::unordered_set<uint64_t> genome{};
    if (expression_by_genome_file != "")
        get_include_set_table(ibf_args, expression_by_genome_file, genome);
    bool const expression_by_genome = (expression_by_genome_file == "");

    // Get expression levels and sizes
    for (unsigned i = 0; i < num_files; i++)
    {
        uint64_t filesize{}; // Store filesize(minimiser_files_given=false) or number of minimisers(minimiser_files_given=true)

        if constexpr(minimiser_files_given)
        {
            uint8_t cutoff;
            read_binary_start(ibf_args, minimiser_files[i], filesize, cutoff);
            cutoffs.push_back(cutoff);
        }
        else
        {
            // Estimate sizes on filesize, assuming every byte translates to one letter (which is obiously not true,
            // because ids contain letters as well), so size might be overestimated. TODO: Find a better estimation!
            unsigned file_iterator = std::accumulate(minimiser_args.samples.begin(), minimiser_args.samples.begin() + i, 0);

            // Determine cutoffs
            if (calculate_cutoffs)
                cutoffs.push_back(calculate_cutoff(minimiser_files[file_iterator], minimiser_args.samples[i]));

            bool const is_compressed = minimiser_files[file_iterator].extension() == ".gz" || minimiser_files[file_iterator].extension() == ".bgzf" || minimiser_files[file_iterator].extension() == ".bz2";
            bool const is_fasta = is_compressed ? check_for_fasta_format(seqan3::format_fasta::file_extensions,minimiser_files[file_iterator].stem())
                                                 : check_for_fasta_format(seqan3::format_fasta::file_extensions, minimiser_files[file_iterator].extension());
            filesize = std::filesystem::file_size(minimiser_files[file_iterator]) * minimiser_args.samples[i] * (is_fasta ? 2 : 1) / (is_compressed ? 1 : 3);
            filesize = filesize/((cutoffs[i] + 1) * (is_fasta ? 1 : 2));
        }
        // If set_expression_thresholds_samplewise is not set the expressions as determined by the first file are used for
        // all files.
        if constexpr (samplewise)
        {
            uint64_t diff{1};
            for (std::size_t c = 0; c < ibf_args.number_expression_thresholds - 1; c++)
            {
                diff = diff * 2;
                sizes[i].push_back(filesize/diff);
            }
            sizes[i].push_back(filesize/diff);
        }
        else if constexpr (minimiser_files_given)
        {
            get_filsize_per_expression_level(minimiser_files[i], ibf_args.number_expression_thresholds, ibf_args.expression_thresholds, sizes[i],
                                             genome, expression_by_genome);
        }
        else
        {
            float diff{1};
            for (std::size_t c = 0; c < ibf_args.number_expression_thresholds - 1; c++)
            {
                diff = ibf_args.expression_thresholds[c+1]/ibf_args.expression_thresholds[c];
                sizes[i].push_back(filesize/diff);
            }
            sizes[i].push_back(filesize/diff);
        }
    }

    std::ofstream outfile_fpr;
    outfile_fpr.open(std::string{ibf_args.path_out} +  "IBF_FPRs.fprs"); // File to store actual false positive rates per experiment.
    // Create IBFs
    std::vector<seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>> ibfs;
    for (unsigned j = 0; j < ibf_args.number_expression_thresholds; j++)
    {
        uint64_t size{0};
        for (unsigned i = 0; i < num_files; i++)
            size = size + sizes[i][j];

        if (size < 1)
        {
            throw std::invalid_argument{std::string("[Error]. The chosen expression threshold is not well picked. If you use the automatic ") +
            std::string("expression threshold determination, please decrease the number of levels. If you use ") +
            std::string("your own expression thresholds, decrease the thresholds from level ") +
            std::to_string(ibf_args.expression_thresholds[j]) +
            std::string(" on.\n")};
        }
        // m = -hn/ln(1-p^(1/h))
        size = static_cast<uint64_t>((-1.0*num_hash*((1.0*size)/num_files))/(std::log(1.0-std::pow(fprs[j], 1.0/num_hash))));
        ibfs.push_back(seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>(
                     seqan3::bin_count{num_files}, seqan3::bin_size{size},
                     seqan3::hash_function_count{num_hash}));

        for (unsigned i = 0; i < num_files; i++)
        {
            double fpr = std::pow(1.0- std::pow(1.0-(1.0/size), num_hash*sizes[i][j]), num_hash);
            outfile_fpr << fpr << " ";
        }
        outfile_fpr << "\n";
    }
    outfile_fpr << "/\n";
    outfile_fpr.close();

    // Add minimisers to ibf
    #pragma omp parallel for schedule(dynamic, chunk_size)
    for (unsigned i = 0; i < num_files; i++)
    {
        robin_hood::unordered_node_map<uint64_t, uint16_t> hash_table{}; // Storage for minimisers
        // Create a smaller cutoff table to save RAM, this cutoff table is only used for constructing the hash table
        // and afterwards discarded.
        robin_hood::unordered_node_map<uint64_t, uint8_t>  cutoff_table;
        std::vector<uint16_t> expression_thresholds;

        // Fill hash table with minimisers.
        if constexpr (minimiser_files_given)
        {
            read_binary(minimiser_files[i], hash_table);
        }
        else
        {
            unsigned file_iterator = std::accumulate(minimiser_args.samples.begin(), minimiser_args.samples.begin() + i, 0);
            for (unsigned f = 0; f < minimiser_args.samples[i]; f++)
            {
               seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{minimiser_files[file_iterator+f]};
               if (minimiser_args.ram_friendly)
                    fill_hash_table_parallel(ibf_args, fin, hash_table, cutoff_table, include_set_table, exclude_set_table,
                               (minimiser_args.include_file != ""), cutoffs[i]);
               else
                    fill_hash_table(ibf_args, fin, hash_table, cutoff_table, include_set_table, exclude_set_table,
                               (minimiser_args.include_file != ""), cutoffs[i]);
            }
            cutoff_table.clear();
        }

        // If set_expression_thresholds_samplewise is not set the expressions as determined by the first file are used for
        // all files.
        if constexpr (samplewise)
        {
           get_expression_thresholds(ibf_args.number_expression_thresholds,
                                 hash_table,
                                 expression_thresholds,
                                 sizes[i],
                                 genome,
                                 cutoffs[i],
                                 expression_by_genome);
           expressions[i] = expression_thresholds;
        }

        // Every minimiser is stored in IBF, if it occurence is greater than or equal to the expression level
        for (auto && elem : hash_table)
        {
            for (int j = ibf_args.number_expression_thresholds - 1; j >= 0 ; --j)
            {
                if constexpr (samplewise)
                {
                    if (elem.second >= expressions[i][j])
                    {
                        ibfs[j].emplace(elem.first, seqan3::bin_index{i});
                        break;
                    }
                }
                else
                {
                    if (elem.second >= ibf_args.expression_thresholds[j])
                    {
                        ibfs[j].emplace(elem.first, seqan3::bin_index{i});
                        break;
                    }
                }
            }
        }
    }

    // Store IBFs
    for (unsigned i = 0; i < ibf_args.number_expression_thresholds; i++)
    {
        std::filesystem::path filename;
        if constexpr(samplewise)
             filename = ibf_args.path_out.string() + "IBF_Level_" + std::to_string(i);
        else
            filename = ibf_args.path_out.string() + "IBF_" + std::to_string(ibf_args.expression_thresholds[i]);

        if (ibf_args.compressed)
        {
            seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf{ibfs[i]};
            store_ibf(ibf, filename);
        }
        else
        {
            store_ibf(ibfs[i], filename);
        }

    }

    // Store all expression thresholds per level.
    if constexpr(samplewise)
    {
        std::ofstream outfile;
        outfile.open(std::string{ibf_args.path_out} +  "IBF_Levels.levels");
        for (unsigned j = 0; j < ibf_args.number_expression_thresholds; j++)
        {
            for (unsigned i = 0; i < num_files; i++)
                 outfile << expressions[i][j] << " ";
            outfile << "\n";
        }
        outfile << "/\n";
        outfile.close();
    }
}

// Create ibfs
std::vector<uint16_t> ibf(std::vector<std::filesystem::path> const & sequence_files,
                          estimate_ibf_arguments & ibf_args, minimiser_arguments & minimiser_args,
                          std::vector<double> & fpr,  std::vector<uint8_t> & cutoffs,
                          std::filesystem::path const expression_by_genome_file, size_t num_hash)
{
    // Declarations
    robin_hood::unordered_node_map<uint64_t, uint16_t> hash_table{}; // Storage for minimisers
    seqan3::concatenated_sequences<seqan3::dna4_vector> sequences; // Storage for sequences in experiment files

    check_cutoffs_samples(sequence_files, minimiser_args.paired, minimiser_args.samples, cutoffs);


    check_expression(ibf_args.expression_thresholds, ibf_args.number_expression_thresholds, expression_by_genome_file);
    check_fpr(ibf_args.number_expression_thresholds, fpr);

    ibf_args.samplewise = (ibf_args.expression_thresholds.size() == 0);

    // Store experiment names
    if (minimiser_args.experiment_names)
    {
        std::ofstream outfile;
        outfile.open(std::string{ibf_args.path_out} + "Stored_Files.txt");
        for (unsigned i = 0; i < minimiser_args.samples.size(); i++)
        {
            outfile  << sequence_files[std::accumulate(minimiser_args.samples.begin(),
                                                minimiser_args.samples.begin()+i, 0)] << "\n";
        }
        outfile.close();
    }

    if (ibf_args.samplewise)
        ibf_helper<true, false>(sequence_files, fpr, ibf_args, cutoffs, num_hash, expression_by_genome_file, minimiser_args);
    else
        ibf_helper<false, false>(sequence_files, fpr, ibf_args, cutoffs, num_hash, expression_by_genome_file, minimiser_args);

    store_args(ibf_args, std::string{ibf_args.path_out} + "IBF_Data");

    return ibf_args.expression_thresholds;
}

// Create ibfs based on the minimiser file
std::vector<uint16_t> ibf(std::vector<std::filesystem::path> const & minimiser_files,
                          estimate_ibf_arguments & ibf_args, std::vector<double> & fpr,
                          std::filesystem::path const expression_by_genome_file,
                          size_t num_hash)
{
    check_expression(ibf_args.expression_thresholds, ibf_args.number_expression_thresholds, expression_by_genome_file);
    check_fpr(ibf_args.number_expression_thresholds, fpr);

    ibf_args.samplewise = (ibf_args.expression_thresholds.size() == 0);

    std::vector<uint8_t> cutoffs{};
    if (ibf_args.samplewise)
        ibf_helper<true>(minimiser_files, fpr, ibf_args, cutoffs, num_hash, expression_by_genome_file);
    else
        ibf_helper<false>(minimiser_files, fpr, ibf_args, cutoffs, num_hash, expression_by_genome_file);

    store_args(ibf_args, std::string{ibf_args.path_out} + "IBF_Data");

    return ibf_args.expression_thresholds;
}

// Actuall minimiser calculation
template<bool parallel = false>
void calculate_minimiser(std::vector<std::filesystem::path> const & sequence_files,
                         robin_hood::unordered_set<uint64_t> const & include_set_table,
                         robin_hood::unordered_set<uint64_t> const & exclude_set_table,
                         min_arguments const & args,
                         minimiser_arguments const & minimiser_args,
                         unsigned const i,
                         std::vector<uint8_t> & cutoffs)
{
    robin_hood::unordered_node_map<uint64_t, uint16_t> hash_table{}; // Storage for minimisers
    uint16_t count{0};
    uint8_t cutoff{0};

    // Create a smaller cutoff table to save RAM, this cutoff table is only used for constructing the hash table
    // and afterwards discarded.
    robin_hood::unordered_node_map<uint64_t, uint8_t>  cutoff_table;
    std::ofstream outfile;
    unsigned file_iterator = std::accumulate(minimiser_args.samples.begin(), minimiser_args.samples.begin() + i, 0);

    bool const calculate_cutoffs = cutoffs.empty();

    if (calculate_cutoffs)
        cutoff = calculate_cutoff(sequence_files[file_iterator], minimiser_args.samples[i]);
    else
        cutoff = cutoffs[i];

    // Fill hash_table with minimisers.
    for (unsigned f = 0; f < minimiser_args.samples[i]; f++)
    {
        seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{sequence_files[file_iterator+f]};
        if constexpr (parallel)
        {
            fill_hash_table_parallel(args, fin, hash_table, cutoff_table, include_set_table, exclude_set_table, (minimiser_args.include_file != ""), cutoff);
        }
        else
        {
            fill_hash_table(args, fin, hash_table, cutoff_table, include_set_table, exclude_set_table, (minimiser_args.include_file != ""), cutoff);
        }
    }
    cutoff_table.clear();

    // Write minimiser and their counts to binary
    outfile.open(std::string{args.path_out} + std::string{sequence_files[file_iterator].stem()}
                 + ".minimiser", std::ios::binary);
    auto hash_size = hash_table.size();

    outfile.write(reinterpret_cast<const char*>(&hash_size), sizeof(hash_size));
    outfile.write(reinterpret_cast<const char*>(&cutoff), sizeof(cutoff));
    outfile.write(reinterpret_cast<const char*>(&args.k), sizeof(args.k));
    outfile.write(reinterpret_cast<const char*>(&args.w_size.get()), sizeof(args.w_size.get()));
    outfile.write(reinterpret_cast<const char*>(&args.s.get()), sizeof(args.s.get()));
    bool ungapped = args.shape.all();
    outfile.write(reinterpret_cast<const char*>(&ungapped), sizeof(ungapped));

    if (!ungapped)
    {
        uint64_t shapesize = args.shape.to_ulong();
        outfile.write(reinterpret_cast<const char*>(&shapesize), sizeof(shapesize));
    }

    for (auto && hash : hash_table)
    {
        outfile.write(reinterpret_cast<const char*>(&hash.first), sizeof(hash.first));
        outfile.write(reinterpret_cast<const char*>(&hash.second), sizeof(hash.second));
    }
    outfile.close();
}

void minimiser(std::vector<std::filesystem::path> const & sequence_files, min_arguments const & args,
               minimiser_arguments & minimiser_args, std::vector<uint8_t> & cutoffs)
{
    // Declarations
    robin_hood::unordered_set<uint64_t> include_set_table{}; // Storage for minimisers in include file
    robin_hood::unordered_set<uint64_t> exclude_set_table{}; // Storage for minimisers in exclude file

    check_cutoffs_samples(sequence_files, minimiser_args.paired, minimiser_args.samples, cutoffs);

    if (minimiser_args.include_file != "")
        get_include_set_table(args, minimiser_args.include_file, include_set_table);
    if (minimiser_args.exclude_file != "")
        get_include_set_table(args, minimiser_args.exclude_file, exclude_set_table);

    seqan3::contrib::bgzf_thread_count = args.threads;
    size_t const chunk_size = std::clamp<size_t>(std::bit_ceil(minimiser_args.samples.size() / args.threads), 1u, 64u);

    // Add minimisers to ibf
    if (minimiser_args.ram_friendly)
    {
        for(unsigned i = 0; i < minimiser_args.samples.size(); i++)
        {
            calculate_minimiser<true>(sequence_files, include_set_table, exclude_set_table, args, minimiser_args, i, cutoffs);
        }
    }
    else
    {
        omp_set_num_threads(args.threads);

        #pragma omp parallel for schedule(dynamic, chunk_size)
        for(unsigned i = 0; i < minimiser_args.samples.size(); i++)
        {
            calculate_minimiser(sequence_files, include_set_table, exclude_set_table, args, minimiser_args, i, cutoffs);
        }
    }

}
