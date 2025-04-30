// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "misc/fill_hash_table.hpp"

// Fill hash table with minimisers greater than the cutoff.
void fill_hash_table(minimiser_arguments const & args,
                     sequence_file_t & fin,
                     robin_hood::unordered_node_map<uint64_t, uint16_t> & hash_table,
                     robin_hood::unordered_node_map<uint64_t, uint8_t> & cutoff_table,
                     robin_hood::unordered_set<uint64_t> const & include_set_table,
                     robin_hood::unordered_set<uint64_t> const & exclude_set_table,
                     bool const only_include,
                     uint8_t cutoff)
{
    // Would result in no minimisers being added to hash_table.
    bool const containment_check_in_empty_include_set = only_include && include_set_table.empty();
    if (containment_check_in_empty_include_set)
        return;

    // Would result in no filtering.
    bool const exclusion_check_in_empty_exclude_set = !only_include && exclude_set_table.empty();

    auto check_value_in_tables = [&](auto && minHash)
    {
        return (only_include && include_set_table.contains(minHash))
            || (!only_include && !exclude_set_table.contains(minHash));
    };

    for (auto & [seq] : fin)
    {
        for (auto && minHash : seqan3::views::minimiser_hash(seq, args.shape, args.w_size, args.s))
        {
            if (exclusion_check_in_empty_exclude_set || check_value_in_tables(minHash))
            {
                // Note that hash_table[minHash] would insert minHash into the map if it does not yet exist.
                auto it = hash_table.find(minHash);

                // If minHash is already in hash table, increase count in hash table.
                if (it != hash_table.end())
                {
                    // Prevent overflow. 65535 is the maximum value for uint16_t.
                    it->second = std::min<uint16_t>(65534u, it->second + 1);
                }
                // If there is no cutoff, add value to hash_table and set count to 1.
                else if (cutoff == 0)
                {
                    hash_table[minHash]++;
                }
                // If minHash equals the cutoff, add it to the hash table and add plus one for the current
                // iteration.
                else if (uint8_t const value = cutoff_table[minHash]; value == cutoff)
                {
                    hash_table[minHash] = value + 1u;
                    cutoff_table.erase(minHash);
                }
                // If none of the above, increase count in cutoff table. Cutoff Table reduces RAM usage by storing
                // minimisers with a low occurence in a smaller hash table.
                else
                {
                    cutoff_table[minHash]++;
                }
            }
        }
    }
}

void fill_hash_table_parallel(minimiser_arguments const & args,
                              sequence_file_t & fin,
                              robin_hood::unordered_node_map<uint64_t, uint16_t> & hash_table,
                              robin_hood::unordered_node_map<uint64_t, uint8_t> & cutoff_table,
                              robin_hood::unordered_set<uint64_t> const & include_set_table,
                              robin_hood::unordered_set<uint64_t> const & exclude_set_table,
                              bool const only_include,
                              uint8_t cutoff)
{
    using value_t = typename robin_hood::unordered_node_map<uint64_t, uint16_t>::value_type;
    using local_hash_table_t = robin_hood::unordered_node_map<uint64_t, std::atomic<uint16_t>>;

    // Step 1: load file in batches
    // Step 2: construct minimiser hash values
    // Step 3: Sort minimiser hash values but pertain their original positions
    // Step 4: Provid simple table with minimiser value and count pair
    // Step 5: Store all non-referenced minimiser values per batch

    size_t thread_count = args.threads;
    // Run thread: block execution and load next 10000 sequences.

    auto seq_file_it = std::ranges::begin(fin);
    using sequence_t = seqan3::dna4_vector;

    auto load_next_chunk = [&]()
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

    auto count_minimiser = [&](auto & local_hash_table, std::vector<uint64_t> minimisers)
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
            auto predicate = [=](uint64_t const other_hash)
            {
                return other_hash == current_minimiser;
            };
            auto next_minimiser_it = std::ranges::find_if_not(minimiser_it, minimiser_end, predicate);
            // minimiser_it now points to the first non equal position
            size_t const minimiser_count = std::ranges::distance(minimiser_it, next_minimiser_it);

            if ((only_include && (include_set_table.contains(current_minimiser)))
                || (!only_include && !exclude_set_table.contains(current_minimiser)))
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

    seqan3::detail::latch sync_point{static_cast<std::ptrdiff_t>(thread_count)};
    seqan3::detail::latch sync_point_2{static_cast<std::ptrdiff_t>(thread_count)};
    std::mutex load_mutex{};
    std::vector<std::pair<size_t, size_t>> intervals{};
    std::optional<seqan3::contrib::fixed_buffer_queue<std::pair<size_t, size_t>>> queue;

    auto job = [&](size_t const thread_id)
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

        { // sequential phase to merge sub tables.
            std::scoped_lock lk{load_mutex};
            if (!is_merged)
            {
                // Merge local hash_tables.
                for (auto & local_hash_table : thread_local_hash_tables)
                {
                    for (auto && [key, counter] : local_hash_table)
                    {
                        if (auto it = hash_table.find(key); it != hash_table.end())
                            it->second =
                                static_cast<uint16_t>(std::min<uint32_t>(65534ul, it->second + counter.load()));
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
                    remaining_minimisers.insert(remaining_minimisers.end(),
                                                local_remaining_minimisers2.begin(),
                                                local_remaining_minimisers2.end());
                }

                std::ranges::sort(remaining_minimisers, std::less{});
                auto minimiser_it = remaining_minimisers.begin();
                auto minimiser_end = remaining_minimisers.end();

                while (minimiser_it != minimiser_end)
                {
                    uint64_t current_minimiser = *minimiser_it;
                    auto predicate = [=](uint64_t const other_hash)
                    {
                        return other_hash == current_minimiser;
                    };
                    auto next_minimiser_it = std::ranges::find_if_not(minimiser_it, minimiser_end, predicate);
                    // minimiser_it now points to the first non equal position
                    size_t const minimiser_count = std::ranges::distance(minimiser_it, next_minimiser_it);

                    if ((only_include && (include_set_table.contains(current_minimiser)))
                        || (!only_include && !exclude_set_table.contains(current_minimiser)))
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
