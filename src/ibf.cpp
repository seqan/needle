#include <chrono>
#include <deque>
#include <iostream>
#include <math.h>
#include <numeric>
#include <string>
#include <algorithm> //reorded because of this error:https://github.com/Homebrew/homebrew-core/issues/44579


#if SEQAN3_WITH_CEREAL
#include <cereal/archives/binary.hpp>
#endif // SEQAN3_WITH_CEREAL

#include <robin_hood.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/io/stream/detail/fast_istreambuf_iterator.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/std/ranges>
#include <seqan3/range/container/concatenated_sequences.hpp>
#include <seqan3/range/container/dynamic_bitset.hpp>
#include <seqan3/range/views/take_until.hpp>

#include "ibf.h"
#include "minimiser.h"

void get_sequences(std::vector<std::filesystem::path> const & sequence_files,
                   seqan3::concatenated_sequences<seqan3::dna4_vector> & sequences, uint16_t min_len, unsigned first,
                   unsigned num_exp)
{
    //Loads all reads from all samples of one experiment to sequences
    //TODO: If not enough memory or too many samples in one experiment, read one file record by record
    for (unsigned i = first; i < (first + num_exp); i++)
    {
        seqan3::sequence_file_input<my_traits> fin{sequence_files[i]};
        for (auto & [seq, id, qual]: fin)
        {
            if (seq.size() >= min_len)
                sequences.push_back(seq);
        }

    }
}

void get_minimisers(arguments const & args, seqan3::concatenated_sequences<seqan3::dna4_vector> const & sequences,
                    robin_hood::unordered_node_map<uint64_t, uint16_t> & hash_table,
                    robin_hood::unordered_set<uint64_t> const & genome_set_table,
                    std::filesystem::path const & include_file, bool only_genome)
{
    // Count minimiser in sequence file
    for (auto seq : sequences)
    {
        for (auto minHash : seqan3::views::minimiser_hash(seq, args.shape, args.w_size, args.s))
        {
            if (!only_genome)
                hash_table[minHash] = std::min<uint16_t>(65534u, hash_table[minHash] + 1);
            //TODO: Use unordered_set contains function instead of find, only works in C++20
            else if ((include_file == "") || (genome_set_table.find(minHash) != genome_set_table.end()))
                hash_table[minHash] = std::min<uint16_t>(65534u, hash_table[minHash] + 1);
        }
    }
}

void fill_hash_table(arguments const & args,
                     seqan3::sequence_file_input<my_traits,  seqan3::fields<seqan3::field::seq>> & fin,
                     robin_hood::unordered_node_map<uint64_t, uint16_t> & hash_table,
                     robin_hood::unordered_set<uint64_t> const & genome_set_table)
{
    for (auto & [seq] : fin)
    {
        for (auto && minHash : seqan3::views::minimiser_hash(seq, args.shape, args.w_size, args.s))
        {
            if (genome_set_table.contains(minHash))
                hash_table[minHash] = std::min<uint16_t>(65534u, hash_table[minHash] + 1);
        }
    }
}

void count(arguments const & args, std::vector<std::filesystem::path> sequence_files, std::filesystem::path genome_file,
           std::filesystem::path out_path, bool paired)
{
    robin_hood::unordered_node_map<uint64_t, uint16_t> hash_table{};
    robin_hood::unordered_set<uint64_t> genome_set_table{};
    std::vector<uint64_t> counter{};
    uint64_t exp{};
    std::ofstream outfile;
    int j;

    seqan3::sequence_file_input<my_traits,  seqan3::fields<seqan3::field::seq>> fin3{genome_file};
    for (auto & [seq] : fin3)
    {
        if (seq.size() >= args.w_size.get())
        {
            for (auto && minHash : seqan3::views::minimiser_hash(seq, args.shape, args.w_size, args.s))
                genome_set_table.insert(minHash);
        }
    }

    seqan3::sequence_file_input<my_traits,  seqan3::fields<seqan3::field::id, seqan3::field::seq>> fin2{genome_file};
    for (unsigned i = 0; i < sequence_files.size(); i++)
    {

        if (paired)
        {
            seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{sequence_files[i]};
            fill_hash_table(args, fin, hash_table, genome_set_table);
            i++;
            fin = sequence_files[i];
            fill_hash_table(args, fin, hash_table, genome_set_table);
        }
        else
        {
            seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{sequence_files[i]};
            fill_hash_table(args, fin, hash_table, genome_set_table);
        }

        outfile.open(std::string{out_path} + std::string{sequence_files[i].stem()} + ".count.out");
        j = 0;
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

// Set arguments that ibf and minimiser use
void set_arguments(ibf_arguments & ibf_args)
{
    // Sort given expression rates
    sort(ibf_args.expression_levels.begin(), ibf_args.expression_levels.end());

    // If the expression levels are not supposed to be set automatically and no number of expression levels are given,
    // throw.
    if (ibf_args.set_expression_levels_samplewise & (ibf_args.number_expression_levels == 0))
    {
        throw std::invalid_argument{"Error. Please set the expression levels or give the number of expression levels."};
    }
    else if (ibf_args.number_expression_levels == 0)
    {
        ibf_args.number_expression_levels = ibf_args.expression_levels.size();
    }
    else if (!ibf_args.set_expression_levels_samplewise & ibf_args.expression_levels.size() == 0)
    {
        std::uint64_t level{2};
        for (std::size_t c = 0; c < ibf_args.number_expression_levels; c++)
        {
            ibf_args.expression_levels.push_back(level);
            level = level *2;
        }
    }
}

void set_include_file(arguments const & args, ibf_arguments & ibf_args,
                   robin_hood::unordered_set<uint64_t> & genome_set_table)
{
    // Generate genome mask
    if (ibf_args.include_file != "")
    {
        seqan3::concatenated_sequences<seqan3::dna4_vector> genome_sequences; // Storage for genome sequences
		get_sequences({ibf_args.include_file}, genome_sequences, args.k);

        // Count minimiser in sequence file
        for (auto seq : genome_sequences)
        {
            for (auto minHash :  seqan3::views::minimiser_hash(seq, args.shape, args.w_size, args.s))
                genome_set_table.insert(minHash);
        }
        genome_sequences.clear();
    }
}

// Set arguments that ibf use
void set_arguments_ibf(arguments const & args, ibf_arguments & ibf_args,
                   robin_hood::unordered_set<uint64_t> & genome_set_table)
{
    if (ibf_args.paired) // If paired is true, a pair is seen as one sample
        ibf_args.samples.assign(ibf_args.sequence_files.size()/2,2);
    if (ibf_args.samples.empty()) // If no samples are given and not paired, every file is seen as one experiment
        ibf_args.samples.assign(ibf_args.sequence_files.size(),1);
    if (ibf_args.cutoffs.empty()) // If no cutoffs are given, every experiment gets a cutoff of zero
        ibf_args.cutoffs.assign(ibf_args.samples.size(),0);
    else if (ibf_args.cutoffs.size() == 1) // If one cutoff is given, every experiment gets this cutoff.
        ibf_args.cutoffs.assign(ibf_args.samples.size(), ibf_args.cutoffs[0]);

    // If sum of ibf_args.samples is not equal to number of files, throw error
    else if (std::accumulate(ibf_args.samples.rbegin(), ibf_args.samples.rend(), 0) != ibf_args.sequence_files.size())
        throw std::invalid_argument{"Error. Incorrect command line input for multiple-samples."};

    set_arguments(ibf_args);

    set_include_file(args, ibf_args, genome_set_table);
}

void check_bin_size(ibf_arguments & ibf_args)
{
    // If no bin size is given or not the right amount, throw error.
    if (ibf_args.bin_size.empty())
        throw std::invalid_argument{"Error. Please give a size for the IBFs in bit."};
    else if (ibf_args.bin_size.size() == 1)
        ibf_args.bin_size.assign(ibf_args.number_expression_levels,ibf_args.bin_size[0]);
    else if (ibf_args.bin_size.size() != ibf_args.number_expression_levels)
        throw std::invalid_argument{"Error. Length of sizes for IBFs in bin_size is not equal to length of expression "
                                    "levels."};
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
void read_header(arguments & args, ibf_arguments & ibf_args, std::filesystem::path filename,
                 std::vector<uint64_t> & counts)
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

    std::ranges::copy(stream_view | seqan3::views::take_until_and_consume(seqan3::is_char<'\n'>),
                                    std::cpp20::back_inserter(buffer));
    ibf_args.cutoffs.push_back(std::stoi(buffer));

    // Read second line = expression levels
    do
    {
        buffer.clear();
        std::ranges::copy(stream_view | seqan3::views::take_until_or_throw(seqan3::is_char<' '> || seqan3::is_char<'\n'>),
                                        std::cpp20::back_inserter(buffer));
        ibf_args.expression_levels.push_back((uint32_t)  std::stoi(buffer));
        if (*stream_it != '\n')
            ++stream_it;
    } while (*stream_it != '\n');
    ++stream_it;

    // Read third line = counts per minimiser per expression level
    do
    {
        buffer.clear();
        std::ranges::copy(stream_view | seqan3::views::take_until_or_throw(seqan3::is_char<' '>),
                                        std::cpp20::back_inserter(buffer));
        counts.push_back(std::stoull(buffer));
        if (*stream_it != '\n')
            ++stream_it;
    } while (*stream_it != '\n');
    fin.close();
}

// Calculate expression levels
void get_expression_levels(arguments const & args, ibf_arguments & ibf_args,
                           robin_hood::unordered_node_map<uint64_t, uint16_t> const & hash_table, unsigned const cutoff,
                           robin_hood::unordered_set<uint64_t> const & genome_set_table)
{
    // Calculate expression levels by taking median recursively
    std::vector<uint16_t> counts;
    for (auto & elem : hash_table)
    {
        if ((elem.second > cutoff) && ((ibf_args.include_file == "")
                                   || (genome_set_table.find(elem.first) != genome_set_table.end())))
            counts.push_back(elem.second);

    }

    // Take the size and get the expression levels by dividing it in equal parts.
    if (cutoff == 0)
    {
        std::size_t prev_pos{0};
        std::size_t num_of_counts = counts.size()/(ibf_args.number_expression_levels + 1);
        for (std::size_t c = 0; c < ibf_args.number_expression_levels; c++)
        {
            std::nth_element(counts.begin() + prev_pos, counts.begin() +  prev_pos + num_of_counts, counts.end());
            prev_pos = prev_pos + num_of_counts;
            ibf_args.expression_levels.push_back(counts[prev_pos]);
        }
    }
    else
    {
        std::size_t prev_pos{0};
        std::size_t num_of_counts = counts.size()/(ibf_args.number_expression_levels);
        // If a cutoff is given, take the start of the counts as first expression
        ibf_args.expression_levels.push_back(*std::min_element(counts.begin(), counts.end()));
        for (std::size_t c = 0; c < ibf_args.number_expression_levels - 1; c++)
        {
            std::nth_element(counts.begin() + prev_pos, counts.begin() +  prev_pos + num_of_counts, counts.end());
            prev_pos = prev_pos + num_of_counts;
            ibf_args.expression_levels.push_back(counts[prev_pos]);
        }
    }
    counts.clear();
}

// Calculates statistics from header files created by minimiser
std::vector<std::tuple<std::vector<uint64_t>, std::vector<uint64_t>>> statistics(arguments & args, ibf_arguments & ibf_args,
std::vector<std::filesystem::path> const & header_files)
{
    // for every expression level a count list of all experiments is created
    std::vector<std::vector<uint64_t>> count_all{};
    std::vector<uint32_t> exp_levels;
    std::vector<std::tuple<std::vector<uint64_t>, std::vector<uint64_t>>> results{};

    // For function call read_header
    ibf_args.expression_levels = {};
    std::vector<uint64_t> counts{};

    uint64_t minimum;
    uint64_t median;
    uint64_t maximum;
    std::vector<uint64_t> first;
    std::vector<uint64_t> second;

    for(auto & file : header_files) // Go over every minimiser file
    {
        read_header(args, ibf_args, file, counts);
        if (count_all.size() == 0) // is true for the very first file
        {
            exp_levels = ibf_args.expression_levels;
            count_all.assign(ibf_args.expression_levels.size(), {});
        }

        for( unsigned i = 0; i < counts.size(); ++i)
            count_all[i].push_back(counts[i]);
        ibf_args.expression_levels.clear();
        counts.clear();
    }

    for( unsigned i = 0; i < exp_levels.size(); ++i)
    {
        minimum = *std::min_element(count_all[i].begin(), count_all[i].end());
        std::nth_element(count_all[i].begin(), count_all[i].begin() + count_all[i].size()/2, count_all[i].end());
        median = count_all[i][count_all[i].size()/2];
        maximum = *std::max_element(count_all[i].begin(), count_all[i].end());

        first = {exp_levels[i]};
        second = {minimum, median, maximum};
        results.push_back(std::make_tuple(first, second));
    }

    return results;

}

std::vector<uint32_t> ibf(arguments const & args, ibf_arguments & ibf_args)
{
    // Declarations
    robin_hood::unordered_node_map<uint64_t, uint16_t> hash_table{}; // Storage for minimisers
    robin_hood::unordered_set<uint64_t> genome_set_table{}; // Storage for minimisers in genome sequences
    seqan3::concatenated_sequences<seqan3::dna4_vector> sequences; // Storage for sequences in experiment files
    std::vector<std::vector<uint32_t>> expressions{};

    set_arguments_ibf(args, ibf_args, genome_set_table);
    check_bin_size(ibf_args);

    if (ibf_args.set_expression_levels_samplewise)
    {
        std::vector<uint32_t> zero_vector(ibf_args.samples.size());
        for (unsigned j = 0; j < ibf_args.number_expression_levels; j++)
            expressions.push_back(zero_vector);
    }

    // Store experiment names
    if (ibf_args.experiment_names)
    {
        std::ofstream outfile;
        outfile.open(std::string{ibf_args.path_out} + "Stored_Files.txt");
        for (unsigned i = 0; i < ibf_args.samples.size(); i++)
        {
            outfile  << ibf_args.sequence_files[std::accumulate(ibf_args.samples.begin(),
                                                ibf_args.samples.begin()+i, 0)] << "\n";
        }
        outfile.close();
    }

    // Create IBFs
    std::vector<seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>> ibfs;
    for (unsigned i = 0; i < ibf_args.number_expression_levels; i++)
        ibfs.push_back(seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>(
                     seqan3::bin_count{ibf_args.samples.size()}, seqan3::bin_size{ibf_args.bin_size[i]},
                     seqan3::hash_function_count{ibf_args.num_hash}));

    // Add minimisers to ibf
    for (unsigned i = 0; i < ibf_args.samples.size(); i++)
    {
        get_sequences(ibf_args.sequence_files, sequences, args.k, std::accumulate(ibf_args.samples.begin(),
                                                                          ibf_args.samples.begin()+i,0),
                                                                          ibf_args.samples[i]);

        get_minimisers(args, sequences, hash_table, genome_set_table, ibf_args.include_file);

        sequences.clear();

        if (ibf_args.set_expression_levels_samplewise)
        {
           get_expression_levels(args,
                                 ibf_args,
                                 hash_table,
                                 ibf_args.cutoffs[i],
                                 genome_set_table);

           for (unsigned k = 0; k < ibf_args.number_expression_levels; k++)
                expressions[k][i] = ibf_args.expression_levels[k];
        }

        // Every minimiser is stored in IBF, if it occurence is greater than or equal to the expression level
        for (auto & elem : hash_table)
        {
            for (int j = ibf_args.number_expression_levels - 1; j >= 0 ; --j)
            {
                if ((ibf_args.expression_levels[j] == 0) & (elem.second > ibf_args.cutoffs[i])) // for comparison with mantis, SBT
                {
                    ibfs[j].emplace(elem.first,seqan3::bin_index{i});
                    break;
                }
                else if (elem.second >= ibf_args.expression_levels[j])
                {
                    ibfs[j].emplace(elem.first,seqan3::bin_index{i});
                    break;
                }
            }
        }
        hash_table.clear();
    }

    // Store IBFs
    for (unsigned i = 0; i < ibf_args.expression_levels.size(); i++)
    {
        std::filesystem::path filename{ibf_args.path_out.string() + "IBF_" + std::to_string(ibf_args.expression_levels[i])};
        if (ibf_args.set_expression_levels_samplewise) // TODO: If this option is choosen the expressions need to be stored
             filename = ibf_args.path_out.string() + "IBF_Level_" + std::to_string(i);
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

    if (ibf_args.set_expression_levels_samplewise)
    {
        std::ofstream outfile;
        outfile.open(std::string{ibf_args.path_out} +  "IBF_Levels.levels");
        for (unsigned j = 0; j < ibf_args.number_expression_levels; j++)
        {
            for (unsigned i = 0; i < ibf_args.samples.size(); i++)
                 outfile << expressions[j][i] << " ";
            outfile << "\n";
        }
        outfile << "/\n";
        outfile.close();
    }

    return ibf_args.expression_levels;
}

// Create ibf based on the minimiser and header files
std::vector<uint32_t> ibf(std::vector<std::filesystem::path> minimiser_files, arguments & args,
                          ibf_arguments & ibf_args)
{
    // Declarations
    robin_hood::unordered_node_map<uint64_t, uint16_t> hash_table{}; // Storage for minimisers
    robin_hood::unordered_set<uint64_t> genome_set_table;
    std::vector<std::vector<uint32_t>> expressions{};

    set_arguments(ibf_args);
    check_bin_size(ibf_args);
    set_include_file(args, ibf_args, genome_set_table);
    if (ibf_args.cutoffs.empty()) // If no cutoffs are given, every experiment gets a cutoff of zero
        ibf_args.cutoffs.assign(minimiser_files.size(),0);
    else if (ibf_args.cutoffs.size() == 1) // If one cutoff is given, every experiment gets this cutoff.
        ibf_args.cutoffs.assign(ibf_args.samples.size(), ibf_args.cutoffs[0]);

    if (ibf_args.set_expression_levels_samplewise)
    {
        std::vector<uint32_t> zero_vector(minimiser_files.size());
        for (unsigned j = 0; j < ibf_args.number_expression_levels; j++)
            expressions.push_back(zero_vector);
    }

    // TODO: if expression levels given does not match the expression levels in header file, get_bin_size gives
    // incorrect results or even an error, if more expression levels are added
    for (unsigned j = 0; j < ibf_args.number_expression_levels; j++)
    {
        seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf =
        seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>(
					   seqan3::bin_count{minimiser_files.size()}, seqan3::bin_size{ibf_args.bin_size[j]},
					   seqan3::hash_function_count{ibf_args.num_hash});

       // Add minimisers to ibf
       for (unsigned i = 0; i < minimiser_files.size(); i++)
       {
           read_binary(hash_table, minimiser_files[i].replace_extension(".minimiser"));

           if (ibf_args.set_expression_levels_samplewise)
           {
              ibf_args.expression_levels.clear();
              get_expression_levels(args,
                                    ibf_args,
                                    hash_table,
                                    ibf_args.cutoffs[i],
                                    genome_set_table);
              for (unsigned k = 0; k < ibf_args.number_expression_levels; k++)
                expressions[k][i] = ibf_args.expression_levels[k];
           }


           // Every minimiser is stored in IBF, if it occurence is greater than or equal to expression level
           for (auto & elem : hash_table)
           {
                 if ((ibf_args.expression_levels[j] == 0) & (elem.second > ibf_args.cutoffs[i])) // for comparison with mantis, SBT
                     ibf.emplace(elem.first,seqan3::bin_index{i});
                 else if ((elem.second >= ibf_args.expression_levels[j]) & (j < ibf_args.number_expression_levels - 1) &(elem.second < ibf_args.expression_levels[j+1]))
                     ibf.emplace(elem.first,seqan3::bin_index{i});
                 else if ((j == ibf_args.number_expression_levels - 1) & (elem.second >= ibf_args.expression_levels[j]))
                     ibf.emplace(elem.first,seqan3::bin_index{i});
           }
           hash_table.clear();
       }

       // Store IBFs
       std::filesystem::path filename{ibf_args.path_out.string() + "IBF_" + std::to_string(ibf_args.expression_levels[j])};
       if (ibf_args.set_expression_levels_samplewise)
           filename = ibf_args.path_out.string() + "IBF_Level_" + std::to_string(j);
       if (args.compressed)
       {
            seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf2{ibf};
		    store_ibf(ibf2, filename);
       }
       else
       {
           store_ibf(ibf, filename);
       }
    }

    if (ibf_args.set_expression_levels_samplewise)
    {
        std::ofstream outfile;
        outfile.open(std::string{ibf_args.path_out} +  "IBF_Levels.levels");
        for (unsigned j = 0; j < ibf_args.number_expression_levels; j++)
        {
            for (unsigned i = 0; i < minimiser_files.size(); i++)
                 outfile << expressions[j][i] << " ";
            outfile << "\n";
        }
        outfile << "/\n";
        outfile.close();
    }

	return ibf_args.expression_levels;
}

void minimiser(arguments const & args, ibf_arguments & ibf_args)
{
    // Declarations
    std::vector<uint32_t> counts;
    robin_hood::unordered_node_map<uint64_t, uint16_t> hash_table{}; // Storage for minimisers
    std::ofstream outfile;
    robin_hood::unordered_set<uint64_t> genome_set_table{}; // Storage for minimisers in genome sequences
    seqan3::concatenated_sequences<seqan3::dna4_vector> sequences; // Storage for sequences in experiment files
    int seen_before{0}; // just to keep track, which sequence files have already been processed, used in stead of:
                        // std::accumulate(ibf_args.samples.begin(), ibf_args.samples.begin()+i,0)

    set_arguments_ibf(args, ibf_args, genome_set_table);

    // Add minimisers to ibf
    for (unsigned i = 0; i < ibf_args.samples.size(); i++)
    {
        get_sequences(ibf_args.sequence_files, sequences, args.k, seen_before, ibf_args.samples[i]);
        get_minimisers(args, sequences, hash_table, genome_set_table, ibf_args.include_file);

        //If no expression values given, determine them
        if (ibf_args.set_expression_levels_samplewise)
        {
           ibf_args.expression_levels.clear();
           get_expression_levels(args,
                                 ibf_args,
                                 hash_table,
                                 ibf_args.cutoffs[i],
                                 genome_set_table);
        }

        counts.assign(ibf_args.expression_levels.size(),0);
        for (auto & elem : hash_table)
        {
            for (int j =  ibf_args.expression_levels.size() - 1; j >= 0; j--)
            {
                if ((ibf_args.expression_levels[j] == 0) & (elem.second > ibf_args.cutoffs[i])) // for comparison with mantis, SBT
                {
                    counts[j]++;
                    break;
                }
                else if ((((elem.second)) >= ibf_args.expression_levels[j]))
                {
                    counts[j]++;
                    break;
                }
            }
        }
        sequences.clear();

        // Write minimiser and their counts to binary
        outfile.open(std::string{ibf_args.path_out} + std::string{ibf_args.sequence_files[seen_before].stem()}
                     + ".minimiser", std::ios::binary);
        for (auto & elem : hash_table)
        {
            outfile.write((char*) &elem.first, sizeof(elem.first));
            outfile.write((char*) &elem.second, sizeof(elem.second));
        }
        outfile.close();
        hash_table.clear();

        // Write header file, containing information about the minimiser counts per expression level
        outfile.open(std::string{ibf_args.path_out} + std::string{ibf_args.sequence_files[seen_before].stem()}
                     + ".header");
        outfile <<  args.s.get() << " " << std::to_string(args.k) << " " << args.w_size.get() << " " << args.shape.to_ulong() << " "
                << ibf_args.cutoffs[i] << "\n";
        for (unsigned k = 0; k < counts.size(); k++)
            outfile  << ibf_args.expression_levels[k] << " ";

        outfile << "\n";
        for (unsigned k = 0; k < counts.size(); k++)
            outfile  << counts[k] << " ";
        outfile << "\n";
        outfile.close();
        seen_before = seen_before + ibf_args.samples[i];
    }

}

void build_ibf(arguments & args, ibf_arguments & ibf_args, float fpr)
{
    std::vector<std::filesystem::path> minimiser_files;
    minimiser(args, ibf_args);
    for (const auto & entry : std::filesystem::directory_iterator(ibf_args.path_out))
    {
        if (entry.path().extension() == ".minimiser")
            minimiser_files.push_back(entry.path());
    }
    // necessary, because std::filesystem::directory_iterator's order is unspecified
    std::sort(minimiser_files.begin(), minimiser_files.end());
    ibf(minimiser_files, args, ibf_args);
    minimiser_files.clear();
}
