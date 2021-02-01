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
#include <seqan3/io/stream/iterator.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/std/ranges>
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
                    robin_hood::unordered_node_map<uint64_t, uint64_t> & hash_table,
                    robin_hood::unordered_set<uint64_t> const & genome_set_table,
                    std::filesystem::path const & include_file, bool only_genome)
{
    // Count minimiser in sequence file
    for (auto seq : sequences)
    {
        for (auto minHash : seqan3::views::minimiser_hash(seq, args.shape, args.w_size, args.s))
        {
            if (!only_genome)
                hash_table[minHash]++;
            //TODO: Use unordered_set contains function instead of find, only works in C++20
            else if ((include_file == "") || (genome_set_table.find(minHash) != genome_set_table.end()))
                hash_table[minHash]++;
        }
    }
}

void count(arguments const & args, std::vector<std::filesystem::path> sequence_files, std::filesystem::path genome_file,
           std::filesystem::path out_path)
{
    seqan3::concatenated_sequences<seqan3::dna4_vector> sequences{};
    robin_hood::unordered_node_map<uint64_t, uint64_t> hash_table{};
    std::vector<std::string> ids{};
    seqan3::concatenated_sequences<seqan3::dna4_vector> genome_sequences{};
    std::vector<uint64_t> counter{};
    uint64_t exp{};
    std::ofstream outfile;
    int j;

    seqan3::sequence_file_input<my_traits> fin2{genome_file};
    for (auto & [seq, id, qual]: fin2)
    {
        if (seq.size() >= args.w_size.get())
        {
            ids.push_back(id);
            genome_sequences.push_back(seq);
        }
    }

    for (unsigned i = 0; i < sequence_files.size(); i++)
    {
        get_sequences(sequence_files, sequences, args.k, i, 1);

        for (auto seq : sequences)
        {
            for (auto minHash : seqan3::views::minimiser_hash(seq, args.shape, args.w_size, args.s))
                hash_table[minHash]++;
        }
        sequences.clear();

        outfile.open(std::string{out_path} + std::string{sequence_files[i].stem()} + ".count.out");
        j = 0;
        for (auto seq : genome_sequences)
        {
            for (auto minHash : seqan3::views::minimiser_hash(seq, args.shape, args.w_size, args.s))
                counter.push_back(hash_table[minHash]);
            std::nth_element(counter.begin(), counter.begin() + counter.size()/2, counter.end());
            exp =  counter[counter.size()/2];
            counter.clear();
            outfile << ids[j] << "\t" << exp << "\n";
            ++j;

        }
        outfile.close();
    }
}

// Set arguments that ibf and minimiser use
void set_arguments(arguments const & args, ibf_arguments & ibf_args,
                   seqan3::concatenated_sequences<seqan3::dna4_vector> & genome_sequences,
                   robin_hood::unordered_set<uint64_t> & genome_set_table)
{
    if (ibf_args.paired) // If paired is true, a pair is seen as one sample
        ibf_args.samples.assign(ibf_args.sequence_files.size()/2,2);
    if (ibf_args.samples.empty()) // If no samples are given and not paired, every file is seen as one experiment
        ibf_args.samples.assign(ibf_args.sequence_files.size(),1);
    if (ibf_args.cutoffs.empty()) // If no cutoffs are given, every experiment gets a cutoff of zero
        ibf_args.cutoffs.assign(ibf_args.samples.size(),0);
    // If sum of ibf_args.samples is not equal to number of files, throw error
    else if (std::accumulate(ibf_args.samples.rbegin(), ibf_args.samples.rend(), 0) != ibf_args.sequence_files.size())
        throw std::invalid_argument{"Error. Incorrect command line input for multiple-samples."};

    // Generate genome mask
    if (ibf_args.include_file != "")
    {
		get_sequences({ibf_args.include_file}, genome_sequences, args.k);

        // Count minimiser in sequence file
        for (auto seq : genome_sequences)
        {
            for (auto minHash :  seqan3::views::minimiser_hash(seq, args.shape, args.w_size, args.s))
                genome_set_table.insert(minHash);
        }
    }

    // Sort given expression rates
    sort(ibf_args.expression_levels.begin(), ibf_args.expression_levels.end());

}

// Reads a binary file minimiser creates
void read_binary(robin_hood::unordered_node_map<uint64_t, uint64_t> & hash_table, std::filesystem::path filename)
{
    std::ifstream fin;
    uint64_t minimiser;
    uint64_t minimiser_count;
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
                 std::vector<uint64_t> & counts, uint32_t & normalized_exp_value)
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
    ibf_args.cutoffs.push_back(std::stoi(buffer));
    buffer.clear();
    std::ranges::copy(stream_view | seqan3::views::take_until_and_consume(seqan3::is_char<' '>),
                                    std::cpp20::back_inserter(buffer));
    ibf_args.normalization_method = buffer;
    buffer.clear();
    std::ranges::copy(stream_view | seqan3::views::take_until_and_consume(seqan3::is_char<'\n'>),
                                    std::cpp20::back_inserter(buffer));
    normalized_exp_value = std::stoi(buffer);

    // Read second line = expression levels
    do
    {
        buffer.clear();
        std::ranges::copy(stream_view | seqan3::views::take_until_or_throw(seqan3::is_char<' '> || seqan3::is_char<'\n'>),
                                        std::cpp20::back_inserter(buffer));
        ibf_args.expression_levels.push_back(std::stof(buffer));
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

// Calculate normalized expression value
uint32_t normalization_method(arguments const & args, ibf_arguments const & ibf_args,
                              seqan3::concatenated_sequences<seqan3::dna4_vector> const & sequences,
                              robin_hood::unordered_node_map<uint64_t, uint64_t> & hash_table, unsigned cutoff,
                              robin_hood::unordered_set<uint64_t> const & genome_set_table)
{
    uint32_t mean; // the normalized expression value

    // Calculate normalized expression value by taking the median of all hash values. Default.
    if (ibf_args.normalization_method == "median")
    {
        std::vector<uint32_t> counts;
        for (auto & elem : hash_table)
        {
            if ((elem.second > cutoff) && ((ibf_args.include_file == "")
                                       || (genome_set_table.find(elem.first) != genome_set_table.end())))
                counts.push_back(elem.second);
        }
        std::nth_element(counts.begin(), counts.begin() + counts.size()/2, counts.end());
        mean = counts[counts.size()/2];
        counts.clear();
    }
    // Calculate normalized expression value by taking the mean of all hash values
    else if (ibf_args.normalization_method == "mean")
    {
        size_t sum_hash{0};
        size_t sum{0};
        for (auto & elem : hash_table)
        {
            if ((elem.second > cutoff) && ((ibf_args.include_file == "")
                                       || (genome_set_table.find(elem.first) != genome_set_table.end())))
            {
                sum = sum + 1;
                sum_hash = sum_hash + elem.second;
            }
        }
        mean = sum_hash/sum;
    }

    return mean;
}

// Calculate expression levels
void get_expression_levels(arguments const & args, ibf_arguments & ibf_args,
                           robin_hood::unordered_node_map<uint64_t, uint64_t> & hash_table, unsigned cutoff,
                           robin_hood::unordered_set<uint64_t> const & genome_set_table, std::size_t filesize)
{
    // Calculate expression levels by taking median recursively
    if (ibf_args.normalization_method == "median")
    {
        std::vector<uint32_t> counts;
        for (auto & elem : hash_table)
        {
            if ((elem.second > cutoff) && ((ibf_args.include_file == "")
                                       || (genome_set_table.find(elem.first) != genome_set_table.end())))
                counts.push_back(elem.second);

        }

        // Start with dev = 4, and prev_pos divided by two, because we ignoare the first median, assuming it is a
        // results of errornous minimizers
        std::size_t dev{4};
        std::size_t prev_pos{counts.size()/2};
        for (std::size_t c = 1; c < ibf_args.number_expression_levels; c++)
        {
            std::nth_element(counts.begin() + prev_pos, counts.begin() +  prev_pos + counts.size()/dev, counts.end());
            prev_pos = prev_pos + counts.size()/dev;
            dev = dev *2;
            ibf_args.expression_levels.push_back(counts[prev_pos]);
        }
        counts.clear();
    }
    // Determine expression levels by file size
    else if (ibf_args.normalization_method == "auto")
    {
        std::array<uint16_t, 4> cutoffs{1, 4, 8, 16};
        std::array<uint64_t, 4> cutoff_bounds{314'572'800, 524'288'000, 1'073'741'824, 3'221'225'472}; // from Mantis/SBT
        std::uint64_t level{32};

        for (size_t k = 0; k < cutoff_bounds.size(); ++k)
        {
            if (filesize <= cutoff_bounds[k])
            {
                level = cutoffs[k];
                break;
            }
        }
        for (std::size_t c = 1; c < ibf_args.number_expression_levels; c++)
        {
            ibf_args.expression_levels.push_back(level);
            level = level *2;
        }
    }
}

// Calculates statistics from header files created by minimiser
std::vector<std::tuple<std::vector<uint64_t>, std::vector<uint64_t>>> statistics(arguments & args, ibf_arguments & ibf_args,
std::vector<std::filesystem::path> const & header_files)
{
    // for every expression level a count list of all experiments is created
    std::vector<std::vector<uint64_t>> count_all{};
    std::vector<uint64_t> normalized_exp_values;
    std::vector<uint64_t> exp_levels;
    std::vector<std::tuple<std::vector<uint64_t>, std::vector<uint64_t>>> results{};

    // For function call read_header
    ibf_args.expression_levels = {};
    std::vector<uint64_t> counts{};
    uint32_t normalized_exp_value{};

    uint64_t avg;
    uint64_t minimum;
    uint64_t median;
    uint64_t maximum;
    std::vector<uint64_t> first;
    std::vector<uint64_t> second;

    for(auto & file : header_files) // Go over every minimiser file
    {
        read_header(args, ibf_args, file, counts, normalized_exp_value);
        if (count_all.size() == 0) // is true for the very first file
        {
            exp_levels = ibf_args.expression_levels;
            count_all.assign(ibf_args.expression_levels.size(), {});
        }

        for( unsigned i = 0; i < counts.size(); ++i)
            count_all[i].push_back(counts[i]);
        normalized_exp_values.push_back(normalized_exp_value);
        ibf_args.expression_levels.clear();
        counts.clear();
    }

    for( unsigned i = 0; i < exp_levels.size(); ++i)
    {
        std::nth_element(normalized_exp_values.begin(), normalized_exp_values.begin() + normalized_exp_values.size()/2,
                         normalized_exp_values.end());
        avg = (1.0 * std::accumulate(normalized_exp_values.begin(),
                                     normalized_exp_values.end(), 0))/normalized_exp_values.size();

        minimum = *std::min_element(count_all[i].begin(), count_all[i].end());
        std::nth_element(count_all[i].begin(), count_all[i].begin() + count_all[i].size()/2, count_all[i].end());
        median = count_all[i][count_all[i].size()/2];
        maximum = *std::max_element(count_all[i].begin(), count_all[i].end());

        first = {exp_levels[i], avg};
        second = {minimum, median, maximum};
        results.push_back(std::make_tuple(first, second));
    }

    return results;

}

std::vector<uint32_t> ibf(arguments const & args, ibf_arguments & ibf_args)
{
    // Declarations
    robin_hood::unordered_node_map<uint64_t, uint64_t> hash_table{}; // Storage for minimisers
    double mean; // the normalized expression value
    std::vector<uint32_t> normal_expression_values;
    seqan3::concatenated_sequences<seqan3::dna4_vector> genome_sequences; // Storage for genome sequences
    robin_hood::unordered_set<uint64_t> genome_set_table{}; // Storage for minimisers in genome sequences
    seqan3::concatenated_sequences<seqan3::dna4_vector> sequences; // Storage for sequences in experiment files

    set_arguments(args, ibf_args, genome_sequences, genome_set_table);

    //If no expression values given, add default
    if (ibf_args.number_expression_levels == 0)
        ibf_args.number_expression_levels = ibf_args.expression_levels.size();

    // If no bin size is given or not the right amount, throw error.
    if (ibf_args.bin_size.empty())
        throw std::invalid_argument{"Error. Please give a size for the IBFs in bit."};
    else if (ibf_args.bin_size.size() == 1)
        ibf_args.bin_size.assign(ibf_args.number_expression_levels,ibf_args.bin_size[0]);
    else if (ibf_args.bin_size.size() != ibf_args.number_expression_levels)
        throw std::invalid_argument{"Error. Length of sizes for IBFs in bin_size is not equal to length of expression "
                                    "levels."};

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
    for (unsigned i = 0; i < ibf_args.expression_levels.size(); i++)
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

        // Calculate normalized expression value in one experiment
        // if genome file is given, calculation are based on genome sequences
        if ((!ibf_args.set_expression_levels) & (ibf_args.include_file != ""))
            mean = normalization_method(args, ibf_args, genome_sequences, hash_table, ibf_args.cutoffs[i], genome_set_table);
        else if ((!ibf_args.set_expression_levels))
            mean = normalization_method(args, ibf_args, sequences, hash_table, ibf_args.cutoffs[i], genome_set_table);
        else
            mean = 1;
        sequences.clear();
        normal_expression_values.push_back(mean);

        if (ibf_args.set_expression_levels)
        {
           get_expression_levels(args,
                                 ibf_args,
                                 hash_table,
                                 ibf_args.cutoffs[i],
                                 genome_set_table,
                                 (ibf_args.samples[i]*std::filesystem::file_size(ibf_args.sequence_files[std::accumulate(ibf_args.samples.begin(),
                                                                                                 ibf_args.samples.begin()+i,0)])));
        }

        // Every minimiser is stored in IBF, if it occurence divided by the mean is greater or equal expression level
        for (int j = ibf_args.number_expression_levels - 1; j >= 0 ; --j)
        {
            for (auto & elem : hash_table)
            {
                if ((ibf_args.expression_levels[j] == 0) & (elem.second > ibf_args.cutoffs[i])) // for comparison with mantis, SBT
                {
                    ibfs[j].emplace(elem.first,seqan3::bin_index{i});
                }
                else if (((elem.second/mean)) >= ibf_args.expression_levels[j])
                {
                    ibfs[j].emplace(elem.first,seqan3::bin_index{i});
                    hash_table.erase(elem.first);
                }
            }
        }
        hash_table.clear();
    }
    genome_sequences.clear();

    // Store IBFs
    for (unsigned i = 0; i < ibf_args.expression_levels.size(); i++)
    {
        std::filesystem::path filename{ibf_args.path_out.string() + "IBF_" + std::to_string(ibf_args.expression_levels[i])};
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

	return normal_expression_values;

}

// Create ibf based on the minimiser and header files
std::vector<uint32_t> ibf(std::vector<std::filesystem::path> minimiser_files, std::filesystem::path header_file,
                          arguments & args, ibf_arguments & ibf_args, float fpr)
{
    // Declarations
    robin_hood::unordered_node_map<uint64_t, uint64_t> hash_table{}; // Storage for minimisers
    uint32_t mean; // the normalized expression value
    std::vector<uint32_t> normal_expression_values;
    // For header file reading
    std::vector<uint64_t> counts{};
    //float normalized_exp_value{};
    std::vector<std::tuple<std::vector<uint64_t>, std::vector<uint64_t>>> statistic_results;
    // Store what user might have entered
    std::vector<uint64_t> expression_levels = ibf_args.expression_levels;
    std::vector<uint64_t> ibf_expression_levels_begin = ibf_args.expression_levels;
    std::string normalization_method = ibf_args.normalization_method;

    if (header_file == "") // If all header files should be considered
    {
        std::vector<std::filesystem::path> header_files{};

        for (auto elem : minimiser_files)
            header_files.push_back(elem.replace_extension(".header"));

        statistic_results = statistics(args, ibf_args, header_files);

        for (unsigned i = 0; i < statistic_results.size(); ++i)
            counts.push_back(std::get<1>(statistic_results[i])[2]); // Get maximum count of all files


        if (expression_levels.size() == 0) // If no expression levels are given
        {
            for (unsigned i = 0; i < statistic_results.size(); ++i)
                expression_levels.push_back(std::get<0>(statistic_results[i])[0]); // Get expression levels
        }
    }
    else // Only given header file is read in
    {
        read_header(args, ibf_args, header_file, counts, mean);

        if (expression_levels.size() == 0) // If no expression levels are given
            expression_levels = ibf_args.expression_levels;
    }

    // TODO: if expression levels given does not match the expression levels in header file, get_bin_size gives
    // incorrect results or even an error, if more expression levels are added
    for (unsigned j = 0; j < expression_levels.size(); j++)
    {
        seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf =
        seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>(
					   seqan3::bin_count{minimiser_files.size()}, seqan3::bin_size{get_bin_size(counts[j], fpr,
                                                                                                ibf_args.num_hash)},
					   seqan3::hash_function_count{ibf_args.num_hash});
       // Add minimisers to ibf
       for (unsigned i = 0; i < minimiser_files.size(); i++)
       {
           read_binary(hash_table, minimiser_files[i].replace_extension(".minimiser"));

           // Get normalized expression value from header file or recalculate it when other method is asked for
           if ((normalization_method == ibf_args.normalization_method) | normalization_method == "")
           {
               read_header(args, ibf_args, minimiser_files[i].replace_extension(".header"), counts, mean);
           }
           // TODO: Add sequences to use - only genome sequences ???
           /*else
           {
               mean = normalization_method(args, ibf_args, sequences, hash_table, ibf_args.cutoffs[i]);
           }
           sequences.clear();*/

           normal_expression_values.push_back(mean);

           // Every minimiser is stored in IBF, if it occurence divided by the mean is greater or equal expression level
           for (auto & elem : hash_table)
           {
                 if ((expression_levels[j] == 0) & (elem.second > ibf_args.cutoffs[i])) // for comparison with mantis, SBT
                     ibf.emplace(elem.first,seqan3::bin_index{i});
                 else if ((((double) elem.second/mean)) >= expression_levels[j])
                     ibf.emplace(elem.first,seqan3::bin_index{i});

           }
           hash_table.clear();
       }
       // Store IBFs
       std::filesystem::path filename{ibf_args.path_out.string() + "IBF_" + std::to_string(expression_levels[j])};
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
    ibf_args.expression_levels = ibf_expression_levels_begin;


	return normal_expression_values;
}

void minimiser(arguments const & args, ibf_arguments & ibf_args)
{
    // Declarations
    std::vector<uint32_t> counts;
    robin_hood::unordered_node_map<uint64_t,uint64_t> hash_table{}; // Storage for minimisers
    std::ofstream outfile;
    uint32_t mean; // the normalized expression value
    seqan3::concatenated_sequences<seqan3::dna4_vector> genome_sequences; // Storage for genome sequences
    robin_hood::unordered_set<uint64_t> genome_set_table{}; // Storage for minimisers in genome sequences
    seqan3::concatenated_sequences<seqan3::dna4_vector> sequences; // Storage for sequences in experiment files
    int seen_before{0}; // just to keep track, which sequence files have already been processed, used in stead of:
                        // std::accumulate(ibf_args.samples.begin(), ibf_args.samples.begin()+i,0)

    set_arguments(args, ibf_args, genome_sequences, genome_set_table);

    // Add minimisers to ibf
    for (unsigned i = 0; i < ibf_args.samples.size(); i++)
    {
        get_sequences(ibf_args.sequence_files, sequences, args.k, seen_before, ibf_args.samples[i]);
        get_minimisers(args, sequences, hash_table, genome_set_table, ibf_args.include_file);

        //If no expression values given, determine them
        if (ibf_args.set_expression_levels)
        {
           get_expression_levels(args,
                                 ibf_args,
                                 hash_table,
                                 ibf_args.cutoffs[i],
                                 genome_set_table,
                                 (ibf_args.samples[i]*std::filesystem::file_size(ibf_args.sequence_files[std::accumulate(ibf_args.samples.begin(),
                                                                                                 ibf_args.samples.begin()+i,0)])));
        }

        // Calculate normalized expression value in one experiment
        // if genome file is given, calculation are based on genome sequences
        if (ibf_args.include_file != "")
            mean = normalization_method(args, ibf_args, genome_sequences, hash_table, ibf_args.cutoffs[i],
                                        genome_set_table);
        else
            mean = normalization_method(args, ibf_args, sequences, hash_table, ibf_args.cutoffs[i], genome_set_table);

        counts.assign(ibf_args.expression_levels.size(),0);
        for (auto & elem : hash_table)
        {
            for (int j =  ibf_args.expression_levels.size() - 1; j >= 0; j--)
            {
                if ((ibf_args.expression_levels[j] == 0) & (elem.second > ibf_args.cutoffs[i])) // for comparison with mantis, SBT
                {
                    counts[j]++;
                }
                else if ((((elem.second/mean)) >= ibf_args.expression_levels[j]))
                {
                    counts[j]++;
                    hash_table.erase(elem.first);
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
                << ibf_args.cutoffs[i] << " " << ibf_args.normalization_method << " " << mean << "\n";
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
    ibf(minimiser_files, "", args, ibf_args, fpr);
    minimiser_files.clear();
}
