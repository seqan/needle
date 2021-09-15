#include <chrono>
#include <deque>
#include <iostream>
#include <math.h>
#include <numeric>
#include <omp.h>
#include <string>
#include <algorithm> //reorded because of this error:https://github.com/Homebrew/homebrew-core/issues/44579

#include <seqan3/std/filesystem>
#include <seqan3/std/ranges>

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
#include <seqan3/utility/container/dynamic_bitset.hpp>

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

uint8_t calculate_cutoff(std::filesystem::path sequence_file, int samples)
{
    // Cutoff according to Mantis paper, divided by two because we store expression levels and
    // -1 because we use "<" and not "<="
    uint16_t const default_cutoff{24};
    uint8_t cutoff{default_cutoff};
    std::array<uint16_t, 4> const cutoffs{0, 1, 4, 9};
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

// Fill hash table with minimisers with cutoff.
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
                // If none of the above, increase count in cutoff table.
                else
                {
                    cutoff_table[minHash]++;
                }
            }
        }
    }
}

void count(min_arguments const & args, std::vector<std::filesystem::path> sequence_files, std::filesystem::path include_file,
           std::filesystem::path exclude_file, bool paired)
{
    robin_hood::unordered_node_map<uint64_t, uint16_t> hash_table{};
    // Create a smaller cutoff table to save RAM, this cutoff table is only used for constructing the hash table
    // and afterwards discarded.
    robin_hood::unordered_node_map<uint64_t, uint8_t>  cutoff_table;
    robin_hood::unordered_set<uint64_t> include_set_table{};
    robin_hood::unordered_set<uint64_t> exclude_set_table{};
    std::vector<uint64_t> counter{};
    uint64_t exp{};
    std::ofstream outfile;
    int j;

    get_include_set_table(args, include_file, include_set_table);
    if (exclude_file != "")
        get_include_set_table(args, exclude_file, exclude_set_table);

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
                 uint64_t & num_of_minimisers)
{
    std::ifstream fin;

    uint32_t window;
    uint64_t buffer;
    fin.open(filename, std::ios::binary);
    fin.read((char*)&buffer, sizeof(buffer));
    num_of_minimisers = buffer;

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
void check_expression(std::vector<uint16_t> & expression_levels, uint8_t & number_expression_levels,
                      std::filesystem::path const expression_by_genome_file)
{
    // Sort given expression rates
    sort(expression_levels.begin(), expression_levels.end());

     // If no expression levels are given and the no number of expression levels is specified, throw.
    if ((number_expression_levels == 0) & (expression_levels.size() == 0))
    {
        throw std::invalid_argument{"Error. Please set the expression levels OR give the number of expression levels."};
    }
    else if ((expression_by_genome_file != "") & (expression_levels.size() > 0))
    {
        throw std::invalid_argument{"Error. The determination of expression levels can not be used with individual levels"
                                    " already given. Please set the expression levels without the option "
                                    "--level-by-genome OR use the number of expression levels with that option."};
    }
    else if (number_expression_levels == 0)
    {
        number_expression_levels = expression_levels.size();
    }
    else if ((number_expression_levels != expression_levels.size()) & (expression_levels.size() > 0))
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

void check_fpr(uint8_t const number_expression_levels, std::vector<double> & fprs)
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

// Calculate expression levels and sizes
void get_expression_levels(uint8_t const number_expression_levels,
                           robin_hood::unordered_node_map<uint64_t, uint16_t> const & hash_table,
                           std::vector<uint16_t> & expression_levels, std::vector<uint64_t> & sizes,
                           robin_hood::unordered_set<uint64_t> const & genome, bool all = true)
{
    // Calculate expression levels by taking median recursively
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
    // First Level
    std::nth_element(counts.begin() + prev_pos, counts.begin() +  prev_pos + counts.size()/dev, counts.end());
    exp = counts[prev_pos + counts.size()/dev];
    prev_pos = prev_pos + counts.size()/dev;
    dev = dev*2;
    expression_levels.push_back(exp);
    sizes.push_back(prev_pos);

    while((expression_levels.size() < number_expression_levels) & (prev_exp < max_elem))
    {
        std::nth_element(counts.begin() + prev_pos, counts.begin() +  prev_pos + counts.size()/dev, counts.end());
        exp = counts[prev_pos + counts.size()/dev];
        prev_pos = prev_pos + counts.size()/dev;
        dev = dev*2;

        if ((exp - prev_exp) > 1)
        {
            expression_levels.push_back(exp);
            sizes.push_back(prev_pos);

        }

        prev_exp = exp;
    }
    while(expression_levels.size() < number_expression_levels)
        expression_levels.push_back(max_elem + 1);
    counts.clear();
}

void get_filsize_per_expression_level(std::filesystem::path filename, uint8_t const number_expression_levels,
                                      std::vector<uint16_t> const & expression_levels, std::vector<uint64_t> & sizes,
                                      robin_hood::unordered_set<uint64_t> const & genome, bool all = true)
{
    std::ifstream fin;
    uint8_t small_buffer;
    uint32_t window;
    uint64_t buffer;
    fin.open(filename, std::ios::binary);
    fin.read((char*)&buffer, sizeof(buffer));
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
    sizes.assign(number_expression_levels, 0);

    while(fin.read((char*)&minimiser, sizeof(minimiser)))
    {
        fin.read((char*)&minimiser_count, sizeof(minimiser_count));
        if (all | genome.contains(minimiser))
        {
            auto p = std::upper_bound(expression_levels.begin(), expression_levels.end(), minimiser_count);
            if(p != expression_levels.begin())
                sizes[(p-expression_levels.begin())-1]++;
            else if (minimiser_count>=expression_levels[number_expression_levels-1])
                sizes[number_expression_levels-1]++;
        }
    }

    fin.close();
}

template<bool samplewise, bool minimiser_files_given = true>
void ibf_helper(std::vector<std::filesystem::path> const & minimiser_files,
                std::vector<double> const & fprs,
                estimate_ibf_arguments & ibf_args, size_t num_hash = 1, std::filesystem::path expression_by_genome_file = "",
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

    bool const calculate_cutoffs = minimiser_args.cutoffs.empty();
    std::vector<uint8_t> file_cutoffs{};

    robin_hood::unordered_set<uint64_t> include_set_table; // Storage for minimisers in include file
    robin_hood::unordered_set<uint64_t> exclude_set_table; // Storage for minimisers in exclude file
    if constexpr(samplewise)
    {
        std::vector<uint16_t> zero_vector(ibf_args.number_expression_levels);
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

    omp_set_num_threads(ibf_args.threads);

    size_t const chunk_size = std::clamp<size_t>(std::bit_ceil(num_files / ibf_args.threads), 8u, 64u);

    // If expression_levels should only be depending on minimsers in a certain genome file, genome is created.
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
            read_binary_start(ibf_args, minimiser_files[i], filesize);
        }
        else
        {
            // Estimate sizes on filesize, assuming every byte translates to one letter (which is obiously not true,
            // because ids contain letters as well), so size might be overestimated
            unsigned file_iterator = std::accumulate(minimiser_args.samples.begin(), minimiser_args.samples.begin() + i, 0);

            // Determine cutoffs
            if (calculate_cutoffs)
                file_cutoffs.push_back(calculate_cutoff(minimiser_files[file_iterator], minimiser_args.samples[i]));

            bool const is_compressed = minimiser_files[file_iterator].extension() == ".gz" || minimiser_files[file_iterator].extension() == ".bgzf" || minimiser_files[file_iterator].extension() == ".bz2";
            bool const is_fasta = is_compressed ? check_for_fasta_format(seqan3::format_fasta::file_extensions,minimiser_files[file_iterator].stem())
                                                 : check_for_fasta_format(seqan3::format_fasta::file_extensions, minimiser_files[file_iterator].extension());
            filesize = std::filesystem::file_size(minimiser_files[file_iterator]) * minimiser_args.samples[i] * (is_fasta ? 2 : 1) / (is_compressed ? 1 : 3);
            if (calculate_cutoffs)
                filesize = filesize/((file_cutoffs[i] + 1) * (is_fasta ? 1 : 2));
            else
                filesize = filesize/((minimiser_args.cutoffs[i] + 1) * (is_fasta ? 1 : 2));
        }
        // If set_expression_levels_samplewise is not set the expressions as determined by the first file are used for
        // all files.
        if constexpr (samplewise)
        {
            uint64_t diff{2};
            for (std::size_t c = 0; c < ibf_args.number_expression_levels - 1; c++)
            {
                diff = diff * 2;
                sizes[i].push_back(filesize/diff);
            }
            sizes[i].push_back(filesize/diff);
        }
        else if constexpr (minimiser_files_given)
        {
            get_filsize_per_expression_level(minimiser_files[i], ibf_args.number_expression_levels, ibf_args.expression_levels, sizes[i],
                                             genome, expression_by_genome);
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

    std::ofstream outfile_fpr;
    outfile_fpr.open(std::string{ibf_args.path_out} +  "IBF_FPRs.fprs");
    // Create IBFs
    std::vector<seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>> ibfs;
    for (unsigned j = 0; j < ibf_args.number_expression_levels; j++)
    {
        uint64_t size{0};
        for (unsigned i = 0; i < num_files; i++)
            size = size + sizes[i][j];

        if (size < 1)
        {
            seqan3::debug_stream << "[Error]. The choosen expression threshold is not well picked. If you use the automatic "
                              << "expression threshold determination, please decrease the number of levels. If you use "
                              << "your own expression thresholds, decrease the thresholds from level "
                              << ibf_args.expression_levels[j] << " on.\n";
        }
        // m = -hn/ln(1-p^(1/h))
        size = static_cast<uint64_t>((-1.0*num_hash*((1.0*size)/num_files))/(std::log(1.0-std::pow(fprs[j], 1.0/num_hash))));
        ibfs.push_back(seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>(
                     seqan3::bin_count{num_files}, seqan3::bin_size{size},
                     seqan3::hash_function_count{num_hash}));

        for (unsigned i = 0; i < num_files; i++)
        {
            double fpr = std::pow(1.0- std::pow(1.0-(1.0/size), num_hash*sizes[i][j]), num_hash);//std::pow((1.0-(std::exp((-1.0*num_hash*sizes[i][j])/size))), num_hash);
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
        std::vector<uint16_t> expression_levels;

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
               if (calculate_cutoffs)
                    fill_hash_table(ibf_args, fin, hash_table, cutoff_table, include_set_table, exclude_set_table,
                                    (minimiser_args.include_file != ""), file_cutoffs[i]);
               else
                    fill_hash_table(ibf_args, fin, hash_table, cutoff_table, include_set_table, exclude_set_table,
                                    (minimiser_args.include_file != ""), minimiser_args.cutoffs[i]);
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
                                 sizes[i],
                                 genome,
                                 expression_by_genome);
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
             filename = ibf_args.path_out.string() + "IBF_Level_" + std::to_string(i);
        else
            filename = ibf_args.path_out.string() + "IBF_" + std::to_string(ibf_args.expression_levels[i]);

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

    if constexpr(samplewise)
    {
        std::ofstream outfile;
        outfile.open(std::string{ibf_args.path_out} +  "IBF_Levels.levels");
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

std::vector<uint16_t> ibf(std::vector<std::filesystem::path> const & sequence_files,
                          estimate_ibf_arguments & ibf_args, minimiser_arguments & minimiser_args,
                          std::vector<double> & fpr,
                          std::filesystem::path const expression_by_genome_file, size_t num_hash)
{
    // Declarations
    robin_hood::unordered_node_map<uint64_t, uint16_t> hash_table{}; // Storage for minimisers
    seqan3::concatenated_sequences<seqan3::dna4_vector> sequences; // Storage for sequences in experiment files

    check_cutoffs_samples(sequence_files, minimiser_args.paired, minimiser_args.samples, minimiser_args.cutoffs);


    check_expression(ibf_args.expression_levels, ibf_args.number_expression_levels, expression_by_genome_file);
    check_fpr(ibf_args.number_expression_levels, fpr);

    ibf_args.samplewise = (ibf_args.expression_levels.size() == 0);

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
        ibf_helper<true, false>(sequence_files, fpr, ibf_args, num_hash, expression_by_genome_file, minimiser_args);
    else
        ibf_helper<false, false>(sequence_files, fpr, ibf_args, num_hash, expression_by_genome_file, minimiser_args);

    store_args(ibf_args, std::string{ibf_args.path_out} + "IBF_Data");

    return ibf_args.expression_levels;
}

// Create ibf based on the minimiser file
std::vector<uint16_t> ibf(std::vector<std::filesystem::path> const & minimiser_files,
                          estimate_ibf_arguments & ibf_args, std::vector<double> & fpr,
                          std::filesystem::path const expression_by_genome_file,
                          size_t num_hash)
{
    check_expression(ibf_args.expression_levels, ibf_args.number_expression_levels, expression_by_genome_file);
    check_fpr(ibf_args.number_expression_levels, fpr);

    ibf_args.samplewise = (ibf_args.expression_levels.size() == 0);

    if (ibf_args.samplewise)
        ibf_helper<true>(minimiser_files, fpr, ibf_args, num_hash, expression_by_genome_file);
    else
        ibf_helper<false>(minimiser_files, fpr, ibf_args, num_hash, expression_by_genome_file);

    store_args(ibf_args, std::string{ibf_args.path_out} + "IBF_Data");

    return ibf_args.expression_levels;
}

void calculate_minimiser(std::vector<std::filesystem::path> const & sequence_files,
                         robin_hood::unordered_set<uint64_t> const & include_set_table,
                         robin_hood::unordered_set<uint64_t> const & exclude_set_table,
                         min_arguments const & args,
                         minimiser_arguments const & minimiser_args,
                         unsigned const i)
{
    robin_hood::unordered_node_map<uint64_t, uint16_t> hash_table{}; // Storage for minimisers
    uint16_t count{0};
    uint8_t cutoff{0};

    // Create a smaller cutoff table to save RAM, this cutoff table is only used for constructing the hash table
    // and afterwards discarded.
    robin_hood::unordered_node_map<uint64_t, uint8_t>  cutoff_table;
    std::ofstream outfile;
    unsigned file_iterator = std::accumulate(minimiser_args.samples.begin(), minimiser_args.samples.begin() + i, 0);

    bool const calculate_cutoffs = minimiser_args.cutoffs.empty();

    if (calculate_cutoffs)
        cutoff = calculate_cutoff(sequence_files[file_iterator], minimiser_args.samples[i]);
    else
        cutoff = minimiser_args.cutoffs[i];

    // Fill hash_table with minimisers.
    for (unsigned f = 0; f < minimiser_args.samples[i]; f++)
    {
        seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{sequence_files[file_iterator+f]};
        fill_hash_table(args, fin, hash_table, cutoff_table, include_set_table, exclude_set_table, (minimiser_args.include_file != ""), cutoff);
    }
    cutoff_table.clear();

    // Write minimiser and their counts to binary
    outfile.open(std::string{args.path_out} + std::string{sequence_files[file_iterator].stem()}
                 + ".minimiser", std::ios::binary);
    auto hash_size = hash_table.size();
    outfile.write(reinterpret_cast<const char*>(&hash_size), sizeof(hash_size));
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

void minimiser(std::vector<std::filesystem::path> const & sequence_files, min_arguments const & args, minimiser_arguments & minimiser_args)
{
    // Declarations
    robin_hood::unordered_set<uint64_t> include_set_table{}; // Storage for minimisers in include file
    robin_hood::unordered_set<uint64_t> exclude_set_table{}; // Storage for minimisers in exclude file

    check_cutoffs_samples(sequence_files, minimiser_args.paired, minimiser_args.samples, minimiser_args.cutoffs);

    if (minimiser_args.include_file != "")
        get_include_set_table(args, minimiser_args.include_file, include_set_table);
    if (minimiser_args.exclude_file != "")
        get_include_set_table(args, minimiser_args.exclude_file, exclude_set_table);

    omp_set_num_threads(args.threads);

    size_t const chunk_size = std::clamp<size_t>(std::bit_ceil(minimiser_args.samples.size() / args.threads),
                                                 1u,
                                                 64u);

    // Add minimisers to ibf
    #pragma omp parallel for schedule(dynamic, chunk_size)
    for(unsigned i = 0; i < minimiser_args.samples.size(); i++)
    {
        calculate_minimiser(sequence_files, include_set_table, exclude_set_table, args, minimiser_args, i);
    }
}
