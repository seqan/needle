#include <chrono>
#include <deque>
#include <iostream>
#include <math.h>
#include <numeric>
#include <omp.h>
#include <string>
#include <algorithm> //reorded because of this error:https://github.com/Homebrew/homebrew-core/issues/44579


#if SEQAN3_WITH_CEREAL
#include <cereal/archives/binary.hpp>
#endif // SEQAN3_WITH_CEREAL

#include <robin_hood.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/algorithm/detail/execution_handler_parallel.hpp>
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
                     robin_hood::unordered_set<uint64_t> const & genome_set_table,
                     bool const only_genome = false, uint8_t cutoff = 0)
{
    // Create a smaller cutoff table to save RAM, this cutoff table is only used for constructing the hash table
    // and afterwards discarded.
    robin_hood::unordered_node_map<uint64_t, uint8_t>  cutoff_table;
    for (auto & [seq] : fin)
    {
        for (auto && minHash : seqan3::views::minimiser_hash(seq, args.shape, args.w_size, args.s))
        {
            if ((only_genome & (genome_set_table.contains(minHash))) | (!only_genome))
            {
                // If minHash is already in hash table, increase count in hash table
                if (hash_table.contains(minHash))
                {
                    hash_table[minHash] = std::min<uint16_t>(65534u, hash_table[minHash] + 1);
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
                    cutoff_table[minHash] = std::min<uint8_t>(254u, cutoff_table[minHash] + 1);
                }
            }
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

    get_include_set_table(args, genome_file, genome_set_table);

    for (unsigned i = 0; i < sequence_files.size(); i++)
    {

        if (paired)
        {
            seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{sequence_files[i]};
            fill_hash_table(args, fin, hash_table, genome_set_table, true);
            i++;
            fin = sequence_files[i];
            fill_hash_table(args, fin, hash_table, genome_set_table, true);
        }
        else
        {
            seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{sequence_files[i]};
            fill_hash_table(args, fin, hash_table, genome_set_table, true);
        }

        outfile.open(std::string{out_path} + std::string{sequence_files[i].stem()} + ".count.out");
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

// Set arguments that ibf and minimiser use
void set_arguments(std::vector<uint16_t> & expression_levels, uint8_t & number_expression_levels,
                   bool const set_expression_levels_samplewise)
{
    // Sort given expression rates
    sort(expression_levels.begin(), expression_levels.end());

     // If no expression levels are given and the no number of expression levels is specified, throw.
    if (((number_expression_levels == 0) & (expression_levels.size() == 0)) |
         (set_expression_levels_samplewise & (expression_levels.size() > 0)))
    {
        throw std::invalid_argument{"Error. Please set the expression levels OR give the number of expression levels."};
    }
    else if (number_expression_levels == 0)
    {
        number_expression_levels = expression_levels.size();
    }
}

// Set arguments that ibf use
void set_arguments_ibf(arguments const & args, std::vector<std::filesystem::path> const & sequence_files,
                       std::filesystem::path const include_file, bool const paired, std::vector<int> & samples,
                       std::vector<uint8_t> & cutoffs, robin_hood::unordered_set<uint64_t> & genome_set_table)
{
    if (paired) // If paired is true, a pair is seen as one sample
        samples.assign(sequence_files.size()/2,2);
    if (samples.empty()) // If no samples are given and not paired, every file is seen as one experiment
        samples.assign(sequence_files.size(),1);
    if (cutoffs.size() == 1) // If one cutoff is given, every experiment gets this cutoff.
        cutoffs.assign(samples.size(), cutoffs[0]);

    // If sum of ibf_args.samples is not equal to number of files, throw error
    else if (std::accumulate(samples.rbegin(), samples.rend(), 0) != sequence_files.size())
        throw std::invalid_argument{"Error. Incorrect command line input for multiple-samples."};

    if (include_file != "")
        get_include_set_table(args, include_file, genome_set_table);
}

void check_bin_size(uint8_t const number_expression_levels, std::vector<size_t> & bin_size)
{
    // If no bin size is given or not the right amount, throw error.
    if (bin_size.empty())
    {
        throw std::invalid_argument{"Error. Please give a size for the IBFs in bit."};
    }
    // If only one ibf size is given, set it for all levels.
    if (bin_size.size() == 1)
    {
        bin_size.assign(number_expression_levels, bin_size[0]);
    }
    else if (bin_size.size() != number_expression_levels)
    {
        throw std::invalid_argument{"Error. Length of sizes for IBFs in bin_size is not equal to length of expression "
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
void read_header(arguments & args, std::vector<uint8_t> & cutoffs, std::filesystem::path filename)
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
    cutoffs.push_back(std::stoi(buffer));

    fin.close();
}

// Calculate expression levels
//TODO: Why do I need args here?
void get_expression_levels(arguments const & args, uint8_t const number_expression_levels,
                           robin_hood::unordered_node_map<uint64_t, uint16_t> const & hash_table,
                           std::vector<uint16_t> & expression_levels)
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
    }
    counts.clear();
}

std::vector<uint16_t> ibf(arguments const & args, ibf_arguments & ibf_args)
{
    // Declarations
    robin_hood::unordered_node_map<uint64_t, uint16_t> hash_table{}; // Storage for minimisers
    robin_hood::unordered_set<uint64_t> genome_set_table{}; // Storage for minimisers in genome sequences
    seqan3::concatenated_sequences<seqan3::dna4_vector> sequences; // Storage for sequences in experiment files
    std::vector<std::vector<uint16_t>> expressions{};

    set_arguments_ibf(args, ibf_args.sequence_files, ibf_args.include_file, ibf_args.paired, ibf_args.samples,
                      ibf_args.cutoffs, genome_set_table);
    if (ibf_args.cutoffs.empty()) // If no cutoffs are given, every experiment gets a cutoff of zero
        ibf_args.cutoffs.assign(ibf_args.samples.size(), 0);

    set_arguments(ibf_args.expression_levels, ibf_args.number_expression_levels,
                  ibf_args.set_expression_levels_samplewise);
    check_bin_size(ibf_args.number_expression_levels, ibf_args.bin_size);

    if (ibf_args.set_expression_levels_samplewise)
    {
        std::vector<uint16_t> zero_vector(ibf_args.samples.size());
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

    unsigned file_iterator{0};
    // Add minimisers to ibf
    for (unsigned i = 0; i < ibf_args.samples.size(); i++)
    {
        // Fill hash table with minimisers.
        for (unsigned f = 0; f < ibf_args.samples[i]; f++)
        {
           seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{ibf_args.sequence_files[file_iterator+f]};
           fill_hash_table(args, fin, hash_table, genome_set_table, (ibf_args.include_file != ""), ibf_args.cutoffs[i]);
        }
        file_iterator = file_iterator + ibf_args.samples[i];

        // If set_expression_levels_samplewise is not set the expressions as determined by the first file are used for
        // all files.
        if (ibf_args.set_expression_levels_samplewise)
        {
           ibf_args.expression_levels.clear();
           get_expression_levels(args,
                                 ibf_args.number_expression_levels,
                                 hash_table,
                                 ibf_args.expression_levels);

           for (unsigned k = 0; k < ibf_args.number_expression_levels; k++)
                expressions[k][i] = ibf_args.expression_levels[k];
        }
        else if (ibf_args.expression_levels.size() == 0)
        {
            get_expression_levels(args,
                                  ibf_args.number_expression_levels,
                                  hash_table,
                                  ibf_args.expression_levels);
        }

        // Every minimiser is stored in IBF, if it occurence is greater than or equal to the expression level
        for (auto && elem : hash_table)
        {
            for (int j = ibf_args.number_expression_levels - 1; j >= 0 ; --j)
            {
                if (elem.second >= ibf_args.expression_levels[j])
                {
                    ibfs[j].emplace(elem.first,seqan3::bin_index{i});
                    break;
                }
            }
        }
        hash_table.clear();
    }

    // Store IBFs
    for (unsigned i = 0; i < ibf_args.number_expression_levels; i++)
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

void fill_ibf(seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> & ibf,
              robin_hood::unordered_node_map<uint64_t, uint16_t> const & hash_table,
              uint32_t const current_expression,
              uint32_t const next_expression,
              unsigned const bin)
{
    // Every minimiser is stored in IBF, if it occurence is greater than or equal to expression level and less
    // than the next expression level.
    for (auto && elem : hash_table)
    {
          if ((elem.second >= current_expression) & (elem.second < next_expression))
              ibf.emplace(elem.first, seqan3::bin_index{bin});
          else if ((next_expression == 0) & (elem.second >= current_expression))
              ibf.emplace(elem.first, seqan3::bin_index{bin});
    }
}

// Create ibf based on the minimiser and header files
//TODO: Why do I need args here?
std::vector<uint16_t> ibf(std::vector<std::filesystem::path> minimiser_files, arguments & args,
                          ibf_arguments & ibf_args)
{

    robin_hood::unordered_node_map<uint64_t, uint16_t> hash_table{}; // Storage for minimisers
    robin_hood::unordered_set<uint64_t> genome_set_table;
    std::vector<std::vector<uint16_t>> expressions{};

    set_arguments(ibf_args.expression_levels, ibf_args.number_expression_levels,
                  ibf_args.set_expression_levels_samplewise);
    check_bin_size(ibf_args.number_expression_levels, ibf_args.bin_size);
    if (ibf_args.include_file != "")
        get_include_set_table(args, ibf_args.include_file, genome_set_table);
    if (ibf_args.cutoffs.empty()) // If no cutoffs are given, every experiment gets a cutoff of zero
        ibf_args.cutoffs.assign(minimiser_files.size(), 0);
    else if (ibf_args.cutoffs.size() == 1) // If one cutoff is given, every experiment gets this cutoff.
        ibf_args.cutoffs.assign(ibf_args.samples.size(), ibf_args.cutoffs[0]);

    if (ibf_args.set_expression_levels_samplewise)
    {
        std::vector<uint16_t> zero_vector(minimiser_files.size());
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
                     seqan3::bin_count{minimiser_files.size()}, seqan3::bin_size{ibf_args.bin_size[i]},
                     seqan3::hash_function_count{ibf_args.num_hash}));

    // Add minimisers to ibf
    for (unsigned i = 0; i < minimiser_files.size(); i++)
    {
        // Fill hash table with minimisers.
        read_binary(hash_table, minimiser_files[i]);

        // If set_expression_levels_samplewise is not set the expressions as determined by the first file are used for
        // all files.
        if (ibf_args.set_expression_levels_samplewise)
        {
           ibf_args.expression_levels.clear();
           get_expression_levels(args,
                                 ibf_args.number_expression_levels,
                                 hash_table,
                                 ibf_args.expression_levels);

           for (unsigned k = 0; k < ibf_args.number_expression_levels; k++)
                expressions[k][i] = ibf_args.expression_levels[k];
        }
        else if (ibf_args.expression_levels.size() == 0)
        {
            get_expression_levels(args,
                                  ibf_args.number_expression_levels,
                                  hash_table,
                                  ibf_args.expression_levels);
        }

        // Every minimiser is stored in IBF, if it occurence is greater than or equal to the expression level
        for (auto && elem : hash_table)
        {
            for (int j = ibf_args.number_expression_levels - 1; j >= 0 ; --j)
            {
                if (elem.second >= ibf_args.expression_levels[j])
                {
                    ibfs[j].emplace(elem.first,seqan3::bin_index{i});
                    break;
                }
            }
        }
        hash_table.clear();
    }

    // Store IBFs
    for (unsigned i = 0; i < ibf_args.number_expression_levels; i++)
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
            for (unsigned i = 0; i < minimiser_files.size(); i++)
                 outfile << expressions[j][i] << " ";
            outfile << "\n";
        }
        outfile << "/\n";
        outfile.close();
    }

    return ibf_args.expression_levels;
}

// Create ibf based on the minimiser and header files
std::vector<uint16_t> ibf_one_ibf_at_a_time(std::vector<std::filesystem::path> minimiser_files, arguments & args,
                          ibf_arguments & ibf_args)
{
    // Declarations
    robin_hood::unordered_node_map<uint64_t, uint16_t> hash_table{}; // Storage for minimisers
    robin_hood::unordered_set<uint64_t> genome_set_table;
    std::vector<std::vector<uint16_t>> expressions{};

    set_arguments(ibf_args.expression_levels, ibf_args.number_expression_levels,
                  ibf_args.set_expression_levels_samplewise);
    check_bin_size(ibf_args.number_expression_levels, ibf_args.bin_size);
    if (ibf_args.include_file != "")
        get_include_set_table(args, ibf_args.include_file, genome_set_table);
    if (ibf_args.cutoffs.empty()) // If no cutoffs are given, every experiment gets a cutoff of zero
        ibf_args.cutoffs.assign(minimiser_files.size(), 0);
    else if (ibf_args.cutoffs.size() == 1) // If one cutoff is given, every experiment gets this cutoff.
        ibf_args.cutoffs.assign(ibf_args.samples.size(), ibf_args.cutoffs[0]);

    if (ibf_args.set_expression_levels_samplewise)
    {
        std::vector<uint16_t> zero_vector(minimiser_files.size());
        for (unsigned j = 0; j < ibf_args.number_expression_levels; j++)
            expressions.push_back(zero_vector);
    }

    // TODO: If expression levels given does not match the expression levels in header file, get_bin_size gives
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
           read_binary(hash_table, minimiser_files[i]);

           if (ibf_args.set_expression_levels_samplewise & (j==0))
           {
              ibf_args.expression_levels.clear();
              get_expression_levels(args,
                                    ibf_args.number_expression_levels,
                                    hash_table,
                                    ibf_args.expression_levels);
              for (unsigned k = 0; k < ibf_args.number_expression_levels; k++)
                 expressions[k][i] = ibf_args.expression_levels[k];
           }
           // If set_expression_levels_samplewise is not set the expressions as determined by the first file are used for
           // all files.
           else if (ibf_args.expression_levels.size() == 0)
           {
               get_expression_levels(args,
                                     ibf_args.number_expression_levels,
                                     hash_table,
                                     ibf_args.expression_levels);
           }

           if (ibf_args.set_expression_levels_samplewise & ((j < ibf_args.number_expression_levels - 1)))
                fill_ibf(ibf, hash_table, expressions[j][i], expressions[j+1][i], i);
           else if (ibf_args.set_expression_levels_samplewise)
                fill_ibf(ibf, hash_table, expressions[j][i], 0, i);
           else if ((j < ibf_args.number_expression_levels - 1))
                fill_ibf(ibf, hash_table, ibf_args.expression_levels[j], ibf_args.expression_levels[j+1], i);
           else
                fill_ibf(ibf, hash_table, ibf_args.expression_levels[j], 0, i);

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

void calculate_minimiser(robin_hood::unordered_set<uint64_t> const & genome_set_table, arguments const & args, ibf_arguments const & ibf_args,
                          unsigned const i)
{
    robin_hood::unordered_node_map<uint64_t, uint16_t> hash_table{}; // Storage for minimisers
    std::ofstream outfile;
    unsigned file_iterator = std::accumulate(ibf_args.samples.begin(), ibf_args.samples.begin() + i, 0);

    bool const calculate_cutoffs = ibf_args.cutoffs.empty();
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
        bool const is_compressed = ibf_args.sequence_files[file_iterator].extension() == ".gz" || ibf_args.sequence_files[file_iterator].extension() == ".bgzf" || ibf_args.sequence_files[file_iterator].extension() == ".bz2";
        bool const is_fasta = is_compressed ? check_for_fasta_format(seqan3::format_fasta::file_extensions, ibf_args.sequence_files[file_iterator].stem())
                                           : check_for_fasta_format(seqan3::format_fasta::file_extensions, ibf_args.sequence_files[file_iterator].extension());
        size_t const filesize = std::filesystem::file_size(ibf_args.sequence_files[file_iterator]) * ibf_args.samples[i] * (is_fasta ? 2 : 1) / (is_compressed ? 1 : 3);

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
        cutoff = ibf_args.cutoffs[i];
    }

    // Fill hash_table with minimisers.
    for (unsigned f = 0; f < ibf_args.samples[i]; f++)
    {
        seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{ibf_args.sequence_files[file_iterator+f]};
        fill_hash_table(args, fin, hash_table, genome_set_table, (ibf_args.include_file != ""), cutoff);
    }

    // Write minimiser and their counts to binary
    outfile.open(std::string{ibf_args.path_out} + std::string{ibf_args.sequence_files[file_iterator].stem()}
                 + ".minimiser", std::ios::binary);

    for (auto && hash : hash_table)
    {
        outfile.write(reinterpret_cast<const char*>(&hash.first), sizeof(hash.first));
        outfile.write(reinterpret_cast<const char*>(&hash.second), sizeof(hash.second));
    }
    outfile.close();
    hash_table.clear();

    // Write header file, containing information about the minimiser counts per expression level
    outfile.open(std::string{ibf_args.path_out} + std::string{ibf_args.sequence_files[file_iterator].stem()}
                 + ".header");
    outfile <<  args.s.get() << " " << std::to_string(args.k) << " " << args.w_size.get() << " " << args.shape.to_ulong() << " "
            << std::to_string(cutoff) << "\n";

    outfile << "\n";
    outfile.close();
}

void minimiser(arguments const & args, ibf_arguments & ibf_args)
{
    // Declarations
    robin_hood::unordered_set<uint64_t> genome_set_table{}; // Storage for minimisers in genome sequences

    set_arguments_ibf(args, ibf_args.sequence_files, ibf_args.include_file, ibf_args.paired, ibf_args.samples,
                      ibf_args.cutoffs, genome_set_table);

    omp_set_num_threads(args.threads);

    // Add minimisers to ibf
    #pragma omp parallel for
    for(unsigned i = 0; i < ibf_args.samples.size(); i++)
    {
        calculate_minimiser(genome_set_table, args, ibf_args, i);
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
