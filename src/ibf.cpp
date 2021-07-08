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

#include <seqan3/std/filesystem>
#include <seqan3/std/ranges>
#include <seqan3/contrib/parallel/buffer_queue.hpp>

#if SEQAN3_WITH_CEREAL
#include <cereal/archives/binary.hpp>
#endif // SEQAN3_WITH_CEREAL

#include <robin_hood.h>

#include <seqan3/alphabet/container/concatenated_sequences.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/io/stream/detail/fast_istreambuf_iterator.hpp>
#include <seqan3/range/views/take_until.hpp>
#include <seqan3/utility/container/dynamic_bitset.hpp>
#include <seqan3/utility/parallel/detail/latch.hpp>

#include "ibf.h"
#include "minimiser.h"

// Create set with hashes from the minimisers from an include file.
void get_include_set_table(arguments const & args, std::filesystem::path const include_file,
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

// Fill hash table with minimisers with cutoff.
void fill_hash_table(arguments const & args,
                     seqan3::sequence_file_input<my_traits,  seqan3::fields<seqan3::field::seq>> & fin,
                     robin_hood::unordered_node_map<uint64_t, uint16_t> & hash_table,
                     robin_hood::unordered_node_map<uint64_t, uint8_t> & cutoff_table,
                     robin_hood::unordered_set<uint64_t> const & genome_set_table,
                     bool const only_genome = false, uint8_t cutoff = 0)
{
    for (auto & [seq] : fin)
    {
        for (auto && minHash : seqan3::views::minimiser_hash(seq, args.shape, args.w_size, args.s))
        {
            if ((only_genome & (genome_set_table.contains(minHash))) | (!only_genome))
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
                // If none of the above, increase count in cutoff table.
                else
                {
                    cutoff_table[minHash]++;
                }
            }
        }
    }
}

void fill_hash_table_parallel(arguments const & args,
                              seqan3::sequence_file_input<my_traits,  seqan3::fields<seqan3::field::seq>> & fin,
                              robin_hood::unordered_node_map<uint64_t, uint16_t> & hash_table,
                              [[maybe_unused]] robin_hood::unordered_node_map<uint64_t, uint8_t> & cutoff_table,
                              [[maybe_unused]] robin_hood::unordered_set<uint64_t> const & genome_set_table,
                              bool const only_genome = false,
                              uint8_t cutoff = 0)
{
    using value_t = typename robin_hood::unordered_node_map<uint64_t, uint16_t>::value_type;
    using local_hash_table_t = robin_hood::unordered_node_map<uint64_t, std::atomic<uint16_t>>;
    // What model do we want to run here?

    // Step 1: load file in batches
    // Step 2: construct minimiser hash values
    // Step 3: Sort minimiser hash values but pertain their original positions
    // Step 4: Provid simple table with minimiser value and count pair
    // Step 5: Store all non-referenced minimiser values per batch

    std::cout << "cutoff = " << (int32_t) cutoff << "\n";
    constexpr size_t thread_count = 1;
    // Run thread: block execution and load next 10000 sequences.

    auto seq_file_it = std::ranges::begin(fin);
    size_t chunk_count{};
    using sequence_t = seqan3::dna4_vector;

    auto load_next_chunk = [&] ()
    {

        constexpr size_t batch_size = 100000;

        std::vector<sequence_t> sequence_batch{};
        sequence_batch.reserve(batch_size);

        // std::cout << "chunk_count = " << chunk_count++ << "\n";
        while (seq_file_it != std::ranges::end(fin) && sequence_batch.size() < batch_size)
        {
            sequence_batch.push_back(seqan3::get<seqan3::field::seq>(*seq_file_it));
            ++seq_file_it;
        }

        return sequence_batch;
    };

    auto count_minimiser = [&] (auto & local_hash_table, std::vector<uint64_t> minimisers)
    {
        std::vector<uint64_t> orphaned_minimiser{};

        // Prepare positional minimiser vector.
        // std::vector<size_t> minimiser_indices{};
        // minimiser_indices.resize(minimisers.size());
        // std::iota(minimiser_indices.begin(), minimiser_indices.end(), 0u);
        // auto projection = [&] (size_t const index) { return minimisers[index]; };

        // std::cout << "Size of minimiser_indices " << minimiser_indices.size() << "\n";
        // seqan3::debug_stream << "minimisers = " << minimisers << "\n";
        // seqan3::debug_stream << "minimiser_indices = " << minimiser_indices << "\n";

        // Sort the minimiser by their value.
        std::ranges::sort(minimisers, std::less{});
        // seqan3::debug_stream << "sorted minimiser_indices = " << minimiser_indices << "\n";

        // Fill the table with all minimiser counts.
        auto minimiser_it = minimisers.begin();
        auto minimiser_end = minimisers.end();

        while (minimiser_it != minimiser_end)
        {
            uint64_t current_minimiser = *minimiser_it;
            auto predicate = [=] (uint64_t const other_hash) { return other_hash == current_minimiser; };
            auto next_minimiser_it = std::ranges::find_if_not(minimiser_it, minimiser_end, predicate);
            // minimiser_it now points to the first non equal position
            size_t const minimiser_count = std::ranges::distance(minimiser_it, next_minimiser_it);

            // if (minimiser_count == 0)
            //     std::cout << '.';
            // If that would be a local hash table, otherwise.
            if ((only_genome & (genome_set_table.contains(current_minimiser))) | (!only_genome))
            {
                if (minimiser_count > cutoff)
                {
                    if (auto it = local_hash_table.find(current_minimiser); it != local_hash_table.end()) // update
                        it->second += minimiser_count; // FIXME: Overflow.
                    else // insert first.
                        local_hash_table[current_minimiser] = minimiser_count;
                }
                else
                { // not above cutoff.
                    // std::cout << minimiser_count << '\t';
                    orphaned_minimiser.insert(orphaned_minimiser.end(), minimiser_it, next_minimiser_it);
                }
            }

            minimiser_it = next_minimiser_it;
        }
        // std::cout << '\n';
        return orphaned_minimiser;
    };

    auto create_queue = [] (size_t const capacity)
    {
        return seqan3::contrib::fixed_buffer_queue<std::pair<size_t, size_t>>{capacity};
    };

    // Block 1:
    constexpr size_t growth_size{100'000'000};
    size_t growth_factor{0};
    std::vector<std::vector<uint64_t>> thread_local_remaining_minimisers{};
    thread_local_remaining_minimisers.resize(thread_count);
    std::vector<local_hash_table_t> thread_local_hash_tables{};
    thread_local_hash_tables.resize(thread_count);
    bool is_merged{false};
    bool create_intervals{true};
    bool fill_queue{true};
    uint64_t min_hash{std::numeric_limits<uint64_t>::max()};
    uint64_t max_hash{0};

    seqan3::detail::latch sync_point{thread_count};
    seqan3::detail::latch sync_point_2{thread_count};
    std::atomic<size_t> remaining_minimisers_size{};
    // std::vector<uint64_t> remaining_minimisers{};
    auto job = [&] (size_t const thread_id)
    {
        std::cout << "start thread " << thread_id << "\n";
        static std::mutex load_mutex{};
        while (true)
        {
            std::vector<sequence_t> sequence_batch{};
            { // critical region
                std::scoped_lock load_lk{load_mutex};
                sequence_batch = load_next_chunk();
            }
            std::cout << "Size of batch " << sequence_batch.size() << "\n";
            if (sequence_batch.empty()) // Stop construction: no more elements are coming
                break;

            // Construct the set of all minimisers for all sequences.
            std::vector<uint64_t> minimisers{};
            minimisers.reserve(sequence_batch.size() * (sequence_batch[0].size() - args.w_size.get() + 1));
            for (auto & sequence : sequence_batch)
                std::ranges::move(sequence | seqan3::views::minimiser_hash(args.shape, args.w_size, args.s),
                                  std::back_inserter(minimisers));

            // std::cout << "Size of minimisers " << minimisers.size() << "\n";
            auto orphaned_minimiser = count_minimiser(thread_local_hash_tables[thread_id], std::move(minimisers));
            std::cout << "Size of orphaned minimisers " << orphaned_minimiser.size() << "\n";
            // remaining_minimiser.reserve(remaining_minimiser.size() + orphaned_minimiser.size());
            thread_local_remaining_minimisers[thread_id].insert(thread_local_remaining_minimisers[thread_id].end(),
                                                                orphaned_minimiser.begin(),
                                                                orphaned_minimiser.end());

            if (thread_local_remaining_minimisers[thread_id].size() > growth_factor * growth_size)
            {
                // Internal reduction step.
                thread_local_remaining_minimisers[thread_id] =
                    count_minimiser(thread_local_hash_tables[thread_id],
                                    std::move(thread_local_remaining_minimisers[thread_id]));
                thread_local_remaining_minimisers[thread_id].reserve((++growth_factor) * growth_size);
            }

            std::cout << "Size of remaining minimisers " << thread_local_remaining_minimisers[thread_id].size() << "\n";
            std::cout << "Size of hash table " << thread_local_hash_tables[thread_id].size() << "\n";
        }

        sync_point.arrive_and_wait();

        // size_t position{};
        // for (uint64_t minimiser : thread_local_remaining_minimisers[thread_id])
        // {
        //     bool is_contained{false};
        //     for (auto & other_hash_tables : thread_local_hash_tables)
        //     {
        //         if (auto it = other_hash_tables.find(minimiser); it != other_hash_tables.end()) // update
        //         {
        //             ++it->second; // FIXME: Overflow.
        //             is_contained = true;
        //             break;
        //         }
        //     }
        //     if (!is_contained)
        //         tmp_remaining_minimisers.push_back(minimiser);

        //     ++position;
        // }

        // std::cout << "old = " <<  thread_local_remaining_minimisers[thread_id].size() << " new = " << tmp_remaining_minimisers.size() << "\n";
        // thread_local_remaining_minimisers[thread_id] = std::move(tmp_remaining_minimisers);

        std::cout << "Here \n";

        {// sequential phase to merge sub tables.
            std::scoped_lock lk{load_mutex};
            std::cout << "Here i\n";
            if (!is_merged)
            {
                std::cout << "Here j\n";
                for (auto & tmp : thread_local_remaining_minimisers)
                {
                    if (tmp.empty())
                        continue;

                    min_hash = std::min(tmp.front(), min_hash);
                    max_hash = std::max(tmp.back(), max_hash);
                }

                std::cout << "Smallest value = " << min_hash << "\n";
                std::cout << "Largest value = " << max_hash << "\n";
                std::cout << "Range value = " << ((max_hash - min_hash + thread_count - 1) / thread_count) << "\n";

                // Merge local hash_tables.
                for (auto & local_hash_table : thread_local_hash_tables)
                {
                    for (auto && [key, counter] : local_hash_table)
                    {
                        if (auto it = hash_table.find(key); it != hash_table.end())
                            it->second += counter.load();
                        else
                            hash_table.insert(value_t{key, counter.load()});
                    }

                    local_hash_table.clear();
                }
                is_merged = true;
            }
        }
        std::cout << "Here 2\n";

        // So we split the range of hash values into number of threads
        if (min_hash > max_hash)
        {
            std::cout << "Early exit\n";
            return;
        }

        constexpr size_t chunk_count = 10;
        size_t const hash_range = ((max_hash - min_hash + (thread_count * chunk_count) - 1) / (thread_count * chunk_count));

        static std::vector<std::pair<size_t, size_t>> intervals{};
        {// sequential phase to merge sub tables.
            std::scoped_lock lk{load_mutex};
            if (create_intervals)
            {
                std::cout << "Preparing the queue\n";
                for (size_t i = 0; i < chunk_count * thread_count; ++i)
                    intervals.emplace_back(min_hash + hash_range * i, min_hash + hash_range * (i + 1));

                create_intervals = false;
                seqan3::debug_stream << "intervals = " << intervals << "\n";
            }
        }
        std::cout << "Here 3\n";

        auto queue = create_queue(intervals.size());
        {// sequential phase to merge sub tables.
            std::scoped_lock lk{load_mutex};
            if (fill_queue)
            {
                for (auto interval : intervals)
                    queue.push(interval);

                queue.close();

                std::cout << "Size of intervals = " << intervals.size() << "\n";
                std::cout << "Size of queue = " << queue.size() << "\n";
                fill_queue = false;
            }
        }
        std::cout << "Here 4\n";
        sync_point_2.arrive_and_wait();

        size_t min_count{};

        // We split the minimiser interval into chunks.
        while (true)
        {
            std::pair<size_t, size_t> interval{};
            if (auto status = queue.wait_pop(interval); status == seqan3::contrib::queue_op_status::closed)
                break;

            auto [thread_range_begin, thread_range_end] = interval;
            std::vector<uint64_t> tmp_remaining_minimisers{};
            for (auto const & tmp : thread_local_remaining_minimisers)
            {
                auto begin_it = std::ranges::lower_bound(tmp, thread_range_begin);
                auto end_it = std::ranges::upper_bound(tmp, thread_range_end);
                min_count += std::ranges::distance(begin_it, end_it);

                tmp_remaining_minimisers.insert(tmp_remaining_minimisers.end(), begin_it, end_it);

                // Now, the remaining minimisers in the set are never reaching the count.
            }
            count_minimiser(thread_local_hash_tables[thread_id], std::move(tmp_remaining_minimisers));
        }

        // sync_point_2.arrive_and_wait();

        {// sequential phase to merge sub tables.
            std::scoped_lock lk{load_mutex};
            // std::cout << "Thread " << thread_id << " [" << thread_range_begin << ", " << thread_range_end << "]\n";
            std::cout << "Counts: " << min_count << "\n";
            std::cout << "close thread " << thread_id << "\n";
        }

    };

    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::thread> thread_pool{};
    for (size_t i = 0; i < thread_count; ++i)
        thread_pool.emplace_back(job, i);

    // Wait for all threads to finish.
    for (auto & thread : thread_pool)
        if (thread.joinable())
            thread.join();

    auto stop = std::chrono::high_resolution_clock::now();
    std::cout << "time initial counting = " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << "s\n";

    // Merge local hash_tables.
    start = std::chrono::high_resolution_clock::now();

    // Record end time
    for (auto & local_hash_table : thread_local_hash_tables)
    {
        for (auto && [key, counter] : local_hash_table)
        {
            if (auto it = hash_table.find(key); it != hash_table.end())
                it->second += counter.load();
            else
                hash_table.insert(value_t{key, counter.load()});
        }

        local_hash_table.clear();
    }
    stop = std::chrono::high_resolution_clock::now();
    std::cout << "time merge ht = " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << "s\n";


    // Reduce minimisers.
    // start = std::chrono::high_resolution_clock::now();
    // std::vector<uint64_t> remaining_minimisers{};
    // for (auto & local_minimisers : thread_local_remaining_minimisers)
    // {
    //     remaining_minimisers.insert(remaining_minimisers.end(), local_minimisers.begin(), local_minimisers.end());
    //     local_minimisers.resize(0);
    // }

    // stop = std::chrono::high_resolution_clock::now();
    // std::cout << "time merge minimisers = " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << "s\n";

    // start = std::chrono::high_resolution_clock::now();
    // auto final_minimiser = count_minimiser(hash_table, std::move(remaining_minimisers));
    // stop = std::chrono::high_resolution_clock::now();
    // std::cout << "time reduce final minimisers = " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << "s\n";
    // std::cout << "Size of final_minimiser " << final_minimiser.size() << "\n";
    std::cout << "Size of hash table " << hash_table.size() << "\n";
    // Final step: lookup all remaining minimiser in the hash table.
    // start = std::chrono::high_resolution_clock::now();
    // for (uint64_t minimiser : final_minimiser)
    // {
    //     if (auto it = hash_table.find(minimiser); it != hash_table.end()) // update
    //         ++it->second; // FIXME: Overflow.
    // }
    // stop = std::chrono::high_resolution_clock::now();
    // std::cout << "time add last minimisers = " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << "s\n";
    // std::cout << "Size of hash table " << hash_table.size() << "\n";
}

void count(arguments const & args, std::vector<std::filesystem::path> sequence_files, std::filesystem::path genome_file,
          bool paired)
{
    robin_hood::unordered_node_map<uint64_t, uint16_t> hash_table{};
    // Create a smaller cutoff table to save RAM, this cutoff table is only used for constructing the hash table
    // and afterwards discarded.
    robin_hood::unordered_node_map<uint64_t, uint8_t>  cutoff_table;
    robin_hood::unordered_set<uint64_t> genome_set_table{};
    std::vector<uint64_t> counter{};
    uint64_t exp{};
    std::ofstream outfile;
    int j;

    get_include_set_table(args, genome_file, genome_set_table);

    for (unsigned i = 0; i < sequence_files.size(); i++)
    {

        if (paired)
        {
            seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{sequence_files[i]};
            fill_hash_table_parallel(args, fin, hash_table, cutoff_table, genome_set_table, true);
            i++;
            fin = sequence_files[i];
            fill_hash_table_parallel(args, fin, hash_table, cutoff_table, genome_set_table, true);
        }
        else
        {
            seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{sequence_files[i]};
            fill_hash_table_parallel(args, fin, hash_table, cutoff_table, genome_set_table, true);
        }
        cutoff_table.clear();

        outfile.open(std::string{args.path_out} + std::string{sequence_files[i].stem()} + ".count.out");
        j = 0;
        seqan3::sequence_file_input<my_traits,  seqan3::fields<seqan3::field::id, seqan3::field::seq>> fin2{genome_file};
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

// Check number of expression levels, sort expression levels
void check_expression(std::vector<uint16_t> & expression_levels, uint8_t & number_expression_levels)
{
    // Sort given expression rates
    sort(expression_levels.begin(), expression_levels.end());

     // If no expression levels are given and the no number of expression levels is specified, throw.
    if ((number_expression_levels == 0) & (expression_levels.size() == 0))
    {
        throw std::invalid_argument{"Error. Please set the expression levels OR give the number of expression levels."};
    }
    else if (number_expression_levels == 0)
    {
        number_expression_levels = expression_levels.size();
    }
}

// Check and set samples and cutoffs
void check_cutoffs_samples(arguments const & args, std::vector<std::filesystem::path> const & sequence_files,
                       std::filesystem::path const include_file, bool const paired, std::vector<int> & samples,
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

void check_fpr(uint8_t const number_expression_levels, std::vector<float> & fprs)
{
    // If no bin size is given or not the right amount, throw error.
    if (fprs.empty())
    {
        throw std::invalid_argument{"Error. Please give a false positive rate for the IBFs."};
    }
    // If only one ibf size is given, set it for all levels.
    if (fprs.size() == 1)
    {
        fprs.assign(number_expression_levels, fprs[0]);
    }
    else if (fprs.size() != number_expression_levels)
    {
        throw std::invalid_argument{"Error. Length of false positive rates for IBFs is not equal to length of expression "
                                    "levels."};
    }
}

// Reads a binary file minimiser creates
void read_binary(robin_hood::unordered_node_map<uint64_t, uint16_t> & hash_table, std::filesystem::path filename)
{
    std::ifstream fin;
    uint64_t minimiser;
    uint16_t minimiser_count;
    fin.open(filename, std::ios::binary);

    while(fin.read((char*)&minimiser, sizeof(minimiser)))
    {
        fin.read((char*)&minimiser_count, sizeof(minimiser_count));
        hash_table[minimiser] = minimiser_count;
    }

    fin.close();
}

// Reads one header file minimiser creates
void read_header(arguments & args, std::vector<uint8_t> & cutoffs, std::filesystem::path filename,
                 uint64_t & num_of_minimisers)
{
    std::ifstream fin;
    fin.open(filename);
    auto stream_view = seqan3::views::istreambuf(fin);
    auto stream_it = std::ranges::begin(stream_view);

    // Read first line
    std::string buffer;
    std::ranges::copy(stream_view | seqan3::views::take_until_and_consume(seqan3::is_char<' '>),
                                    std::cpp20::back_inserter(buffer));
    args.s = seqan3::seed{(uint64_t) std::stoull(buffer)};
    buffer.clear();
    std::ranges::copy(stream_view | seqan3::views::take_until_and_consume(seqan3::is_char<' '>),
                                    std::cpp20::back_inserter(buffer));
    args.k = std::stoi(buffer);
    buffer.clear();
    std::ranges::copy(stream_view | seqan3::views::take_until_and_consume(seqan3::is_char<' '>),
                                    std::cpp20::back_inserter(buffer));
    args.w_size = seqan3::window_size{(uint32_t) std::stoi(buffer)};
    buffer.clear();
    std::ranges::copy(stream_view | seqan3::views::take_until_and_consume(seqan3::is_char<' '>),
                                    std::cpp20::back_inserter(buffer));
    uint64_t shape = (uint64_t) std::stoull(buffer);
    buffer.clear();
    if (shape == 0)
            args.shape = seqan3::ungapped{args.k};
        else
            args.shape = seqan3::bin_literal{shape};

    std::ranges::copy(stream_view | seqan3::views::take_until_and_consume(seqan3::is_char<' '>),
                                    std::cpp20::back_inserter(buffer));
    cutoffs.push_back(std::stoi(buffer));
    buffer.clear();

    std::ranges::copy(stream_view | seqan3::views::take_until_and_consume(seqan3::is_char<'\n'>),
                                    std::cpp20::back_inserter(buffer));
    num_of_minimisers = (uint64_t) std::stoull(buffer);

    fin.close();
}

// Calculate expression levels and sizes
void get_expression_levels(uint8_t const number_expression_levels,
                           robin_hood::unordered_node_map<uint64_t, uint16_t> const & hash_table,
                           std::vector<uint16_t> & expression_levels, std::vector<uint64_t> & sizes)
{
    // Calculate expression levels by taking median recursively
    std::vector<uint16_t> counts;
    for (auto && elem : hash_table)
    {
        counts.push_back(elem.second);
    }

    std::size_t dev{2};
    std::size_t prev_pos{0};
    for (std::size_t c = 0; c < number_expression_levels; c++)
    {
        std::nth_element(counts.begin() + prev_pos, counts.begin() +  prev_pos + counts.size()/dev, counts.end());
        prev_pos = prev_pos + counts.size()/dev;
        dev = dev *2;
        expression_levels.push_back(counts[prev_pos]);
        sizes.push_back(prev_pos);
    }
    counts.clear();
}

void get_sizes(uint8_t const number_expression_levels,
               robin_hood::unordered_node_map<uint64_t, uint16_t> const & hash_table,
               std::vector<uint16_t> const & expression_levels, std::vector<uint64_t> & sizes)
{
    std::vector<uint16_t> counts;
    for (auto && elem : hash_table)
    {
        counts.push_back(elem.second);
    }

    uint64_t prev_size{0};
    for (std::size_t c = 0; c < number_expression_levels; c++)
    {
        auto it = std::find(counts.begin(), counts.end(), expression_levels[c]);
        uint64_t index = static_cast<uint64_t>(*it);
        prev_size = counts.size() - std::count(counts.begin(), counts.end(), expression_levels[c]) + index - prev_size;
        sizes.push_back(prev_size);
    }
}

template<bool samplewise, bool minimiser_files_given = true>
void ibf_helper(std::vector<std::filesystem::path> const & minimiser_files, arguments const & args,
                ibf_arguments const & ibf_args, minimiser_arguments const & minimiser_args = {})
{
    size_t num_files;
    if constexpr (minimiser_files_given)
        num_files = minimiser_files.size();
    else
        num_files = minimiser_args.samples.size();

    std::vector<std::vector<uint16_t>> expressions{};
    std::vector<std::vector<uint64_t>> sizes{};
    sizes.assign(num_files, {});

    robin_hood::unordered_set<uint64_t> genome_set_table; // Storage for minimisers in genome sequences
    if constexpr(samplewise)
    {
        std::vector<uint16_t> zero_vector(ibf_args.number_expression_levels);
        for (unsigned j = 0; j < num_files; j++)
            expressions.push_back(zero_vector);

    }
    if constexpr (!minimiser_files_given)
    {
        if (minimiser_args.include_file != "")
            get_include_set_table(args, minimiser_args.include_file, genome_set_table);
    }

    omp_set_num_threads(args.threads);

    size_t const chunk_size = std::clamp<size_t>(std::bit_ceil(num_files / args.threads),
                                                 8u,
                                                 64u);
    // Get expression levels and sizes
    #pragma omp parallel for schedule(dynamic, chunk_size)
    for (unsigned i = 0; i < num_files; i++)
    {
        uint64_t filesize{}; // Store filesize(minimiser_files_given=false) or number of minimisers(minimiser_files_given=true)

        if constexpr(minimiser_files_given)
        {
            arguments args2{};
            std::vector<uint8_t> cutoffs2{};
            read_header(args2, cutoffs2, std::string{minimiser_files[i].parent_path()} + "/" + std::string{minimiser_files[i].stem()} + ".header", filesize);
        }
        else
        {
            // Estimate sizes on filesize, assuming every byte translates to one letter (which is obiously not true,
            // because ids contain letters as well), so size might be overestimated
            unsigned file_iterator = std::accumulate(minimiser_args.samples.begin(), minimiser_args.samples.begin() + i, 0);
            bool const is_compressed = minimiser_files[file_iterator].extension() == ".gz" || minimiser_files[file_iterator].extension() == ".bgzf" || minimiser_files[file_iterator].extension() == ".bz2";
            bool const is_fasta = is_compressed ? check_for_fasta_format(seqan3::format_fasta::file_extensions,minimiser_files[file_iterator].stem())
                                                 : check_for_fasta_format(seqan3::format_fasta::file_extensions, minimiser_files[file_iterator].extension());
            filesize = std::filesystem::file_size(minimiser_files[file_iterator]) * minimiser_args.samples[i] * (is_fasta ? 2 : 1) / (is_compressed ? 1 : 3);
        }
        // If set_expression_levels_samplewise is not set the expressions as determined by the first file are used for
        // all files.
        if constexpr (samplewise)
        {
            uint64_t diff{1};
            for (std::size_t c = 0; c < ibf_args.number_expression_levels - 1; c++)
            {
                diff = diff * 2;
                sizes[i].push_back(filesize/diff);
            }
            sizes[i].push_back(filesize/diff);
        }
        else
        {
            float diff{1};
            for (std::size_t c = 0; c < ibf_args.number_expression_levels - 1; c++)
            {
                diff = ibf_args.expression_levels[c+1]/ibf_args.expression_levels[c];
                sizes[i].push_back(filesize/diff);
            }
            sizes[i].push_back(filesize/diff);

        }
    }

    // Create IBFs
    std::vector<seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>> ibfs;
    for (unsigned j = 0; j < ibf_args.number_expression_levels; j++)
    {
        uint64_t size{0};
        for (unsigned i = 0; i < num_files; i++)
            size = size + sizes[i][j];
        // m = -hn/ln(1-p^(1/h))
        size = static_cast<uint64_t>((-1.0*ibf_args.num_hash*((1.0*size)/num_files))/(std::log(1.0-std::pow(ibf_args.fpr[j], 1.0/ibf_args.num_hash))));
        ibfs.push_back(seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>(
                     seqan3::bin_count{num_files}, seqan3::bin_size{size},
                     seqan3::hash_function_count{ibf_args.num_hash}));
    }

    // Add minimisers to ibf
    #pragma omp parallel for schedule(dynamic, chunk_size)
    for (unsigned i = 0; i < num_files; i++)
    {
        robin_hood::unordered_node_map<uint64_t, uint16_t> hash_table{}; // Storage for minimisers
        robin_hood::unordered_node_map<uint64_t, uint16_t> hash_table_test{}; // Storage for minimisers
        // Create a smaller cutoff table to save RAM, this cutoff table is only used for constructing the hash table
        // and afterwards discarded.
        robin_hood::unordered_node_map<uint64_t, uint8_t>  cutoff_table;
        std::vector<uint16_t> expression_levels;

        // Fill hash table with minimisers.
        if constexpr (minimiser_files_given)
        {
            read_binary(hash_table, minimiser_files[i]);
        }
        else
        {
            unsigned file_iterator = std::accumulate(minimiser_args.samples.begin(), minimiser_args.samples.begin() + i, 0);
            for (unsigned f = 0; f < minimiser_args.samples[i]; f++)
            {
               seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{minimiser_files[file_iterator+f]};
               fill_hash_table_parallel(args, fin, hash_table, cutoff_table, genome_set_table, (minimiser_args.include_file != ""), minimiser_args.cutoffs[i]);
            }
            cutoff_table.clear();
        }

        // If set_expression_levels_samplewise is not set the expressions as determined by the first file are used for
        // all files.
        if constexpr (samplewise)
        {
           get_expression_levels(ibf_args.number_expression_levels,
                                 hash_table,
                                 expression_levels,
                                 sizes[i]);
           expressions[i] = expression_levels;
        }

        // Every minimiser is stored in IBF, if it occurence is greater than or equal to the expression level
        for (auto && elem : hash_table)
        {
            for (int j = ibf_args.number_expression_levels - 1; j >= 0 ; --j)
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
                    if (elem.second >= ibf_args.expression_levels[j])
                    {
                        ibfs[j].emplace(elem.first, seqan3::bin_index{i});
                        break;
                    }
                }
            }
        }
    }

    // Store IBFs
    for (unsigned i = 0; i < ibf_args.number_expression_levels; i++)
    {
        std::filesystem::path filename;
        if constexpr(samplewise)
             filename = args.path_out.string() + "IBF_Level_" + std::to_string(i);
        else
            filename = args.path_out.string() + "IBF_" + std::to_string(ibf_args.expression_levels[i]);

        if (args.compressed)
        {
            seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf{ibfs[i]};
            store_ibf(ibf, filename);
        }
        else
        {
            store_ibf(ibfs[i], filename);
        }

    }

    if constexpr(samplewise)
    {
        std::ofstream outfile;
        outfile.open(std::string{args.path_out} +  "IBF_Levels.levels");
        for (unsigned j = 0; j < ibf_args.number_expression_levels; j++)
        {
            for (unsigned i = 0; i < num_files; i++)
                 outfile << expressions[i][j] << " ";
            outfile << "\n";
        }
        outfile << "/\n";
        outfile.close();
    }
}

std::vector<uint16_t> ibf(std::vector<std::filesystem::path> const & sequence_files, arguments const & args,
                          ibf_arguments & ibf_args, minimiser_arguments & minimiser_args)
{
    // Declarations
    robin_hood::unordered_node_map<uint64_t, uint16_t> hash_table{}; // Storage for minimisers
    seqan3::concatenated_sequences<seqan3::dna4_vector> sequences; // Storage for sequences in experiment files

    check_cutoffs_samples(args, sequence_files, minimiser_args.include_file, minimiser_args.paired,
                      minimiser_args.samples, minimiser_args.cutoffs);
    if (minimiser_args.cutoffs.empty()) // If no cutoffs are given, every experiment gets a cutoff of zero
        minimiser_args.cutoffs.assign(minimiser_args.samples.size(), 0);

    check_expression(ibf_args.expression_levels, ibf_args.number_expression_levels);
    check_fpr(ibf_args.number_expression_levels, ibf_args.fpr);

    bool samplewise = (ibf_args.expression_levels.size() == 0);

    // Store experiment names
    if (minimiser_args.experiment_names)
    {
        std::ofstream outfile;
        outfile.open(std::string{args.path_out} + "Stored_Files.txt");
        for (unsigned i = 0; i < minimiser_args.samples.size(); i++)
        {
            outfile  << sequence_files[std::accumulate(minimiser_args.samples.begin(),
                                                minimiser_args.samples.begin()+i, 0)] << "\n";
        }
        outfile.close();
    }

    if (samplewise)
        ibf_helper<true, false>(sequence_files, args, ibf_args, minimiser_args);
    else
        ibf_helper<false, false>(sequence_files, args, ibf_args, minimiser_args);

    return ibf_args.expression_levels;
}

// Create ibf based on the minimiser and header files
std::vector<uint16_t> ibf(std::vector<std::filesystem::path> const & minimiser_files, arguments const & args,
                          ibf_arguments & ibf_args)
{
    check_expression(ibf_args.expression_levels, ibf_args.number_expression_levels);
    check_fpr(ibf_args.number_expression_levels, ibf_args.fpr);

    bool const samplewise = (ibf_args.expression_levels.size() == 0);

    if (samplewise)
        ibf_helper<true>(minimiser_files, args, ibf_args);
    else
        ibf_helper<false>(minimiser_files, args, ibf_args);

    return ibf_args.expression_levels;
}

void calculate_minimiser(std::vector<std::filesystem::path> const & sequence_files,
                         robin_hood::unordered_set<uint64_t> const & genome_set_table,
                         arguments const & args,
                         minimiser_arguments const & minimiser_args,
                         unsigned const i)
{
    robin_hood::unordered_node_map<uint64_t, uint16_t> hash_table{}; // Storage for minimisers
    // Create a smaller cutoff table to save RAM, this cutoff table is only used for constructing the hash table
    // and afterwards discarded.
    robin_hood::unordered_node_map<uint64_t, uint8_t>  cutoff_table;
    std::ofstream outfile;
    unsigned file_iterator = std::accumulate(minimiser_args.samples.begin(), minimiser_args.samples.begin() + i, 0);

    bool const calculate_cutoffs = minimiser_args.cutoffs.empty();
    // Cutoff according to Mantis paper
    uint16_t const default_cutoff{50};
    std::array<uint16_t, 4> const cutoffs{1, 3, 10, 20};
    std::array<uint64_t, 4> const cutoff_bounds{314'572'800, 524'288'000, 1'073'741'824, 3'221'225'472};
    uint16_t cutoff{default_cutoff};

    if (calculate_cutoffs)
    {
        uint16_t count{0};
        // Since the curoffs are based on the filesize of a gzipped fastq file, we try account for the other cases:
        // We multiply by two if we have fasta input.
        // We divide by 3 if the input is not compressed.
        bool const is_compressed = sequence_files[file_iterator].extension() == ".gz" || sequence_files[file_iterator].extension() == ".bgzf" || sequence_files[file_iterator].extension() == ".bz2";
        bool const is_fasta = is_compressed ? check_for_fasta_format(seqan3::format_fasta::file_extensions, sequence_files[file_iterator].stem())
                                           : check_for_fasta_format(seqan3::format_fasta::file_extensions, sequence_files[file_iterator].extension());
        size_t const filesize = std::filesystem::file_size(sequence_files[file_iterator]) * minimiser_args.samples[i] * (is_fasta ? 2 : 1) / (is_compressed ? 1 : 3);

        for (size_t k = 0; k < cutoff_bounds.size(); ++k)
        {
            if (filesize <= cutoff_bounds[k])
            {
                cutoff = cutoffs[k];
                break;
            }
        }
    }
    else
    {
        cutoff = minimiser_args.cutoffs[i];
    }

    // Fill hash_table with minimisers.
    for (unsigned f = 0; f < minimiser_args.samples[i]; f++)
    {
        seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{sequence_files[file_iterator+f]};
        fill_hash_table_parallel(args, fin, hash_table, cutoff_table, genome_set_table, (minimiser_args.include_file != ""), cutoff);
    }
    cutoff_table.clear();

    // Write minimiser and their counts to binary
    outfile.open(std::string{args.path_out} + std::string{sequence_files[file_iterator].stem()}
                 + ".minimiser", std::ios::binary);

    for (auto && hash : hash_table)
    {
        outfile.write(reinterpret_cast<const char*>(&hash.first), sizeof(hash.first));
        outfile.write(reinterpret_cast<const char*>(&hash.second), sizeof(hash.second));
    }
    outfile.close();

    // Write header file, containing information about the minimiser counts per expression level
    outfile.open(std::string{args.path_out} + std::string{sequence_files[file_iterator].stem()}
                 + ".header");
    outfile <<  args.s.get() << " " << std::to_string(args.k) << " " << args.w_size.get() << " " << args.shape.to_ulong() << " "
            << std::to_string(cutoff) << " " << std::to_string(hash_table.size()) <<"\n";

    outfile << "\n";
    hash_table.clear();
    outfile.close();
}

void minimiser(std::vector<std::filesystem::path> const & sequence_files, arguments const & args, minimiser_arguments & minimiser_args)
{
    // Declarations
    robin_hood::unordered_set<uint64_t> genome_set_table{}; // Storage for minimisers in genome sequences

    check_cutoffs_samples(args, sequence_files, minimiser_args.include_file, minimiser_args.paired, minimiser_args.samples,
                      minimiser_args.cutoffs);

    if (minimiser_args.include_file != "")
        get_include_set_table(args, minimiser_args.include_file, genome_set_table);

    omp_set_num_threads(args.threads);

    size_t const chunk_size = std::clamp<size_t>(std::bit_ceil(minimiser_args.samples.size() / args.threads),
                                                 1u,
                                                 64u);

    // Add minimisers to ibf
    #pragma omp parallel for schedule(dynamic, chunk_size)
    for(unsigned i = 0; i < minimiser_args.samples.size(); i++)
    {
        calculate_minimiser(sequence_files, genome_set_table, args, minimiser_args, i);
    }
}
