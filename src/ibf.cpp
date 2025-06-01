// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "ibf.hpp"

#include <omp.h>

#include <cereal/archives/binary.hpp>

#include <hibf/build/bin_size_in_bits.hpp>

#include "misc/calculate_cutoff.hpp"
#include "misc/check_cutoffs_samples.hpp"
#include "misc/check_for_fasta_format.hpp"
#include "misc/filenames.hpp"
#include "misc/fill_hash_table.hpp"
#include "misc/get_expression_thresholds.hpp"
#include "misc/get_include_set_table.hpp"
#include "misc/stream.hpp"

// Check number of expression levels, sort expression levels
void check_expression(std::vector<uint16_t> & expression_thresholds,
                      uint8_t & number_expression_thresholds,
                      std::filesystem::path const & expression_by_genome_file)
{
    // Sort given expression rates
    std::ranges::sort(expression_thresholds);

    // If no expression levels are given and the no number of expression levels is specified, throw.
    if ((number_expression_thresholds == 0) && (expression_thresholds.size() == 0))
    {
        throw std::invalid_argument{"Please set the expression levels OR give the number of expression levels."};
    }
    else if ((expression_by_genome_file != "") && (expression_thresholds.size() > 0))
    {
        throw std::invalid_argument{"The determination of expression levels can not be used with individual levels"
                                    " already given. Please set the expression levels without the option "
                                    "--level-by-genome OR use the number of expression levels with that option."};
    }
    else if (number_expression_thresholds == 0)
    {
        number_expression_thresholds = expression_thresholds.size();
    }
    else if ((number_expression_thresholds != expression_thresholds.size()) && (expression_thresholds.size() > 0))
    {
        throw std::invalid_argument{"Please set the expression levels OR give the number of expression levels."};
    }
}

// Check input of fpr
void check_fpr(uint8_t const number_expression_thresholds, std::vector<double> & fprs)
{
    // If no bin size is given or not the right amount, throw error.
    if (fprs.empty())
    {
        throw std::invalid_argument{"Please give a false positive rate for the IBFs."};
    }
    // If only one ibf size is given, set it for all thresholds.
    if (fprs.size() == 1)
    {
        double const fpr = fprs[0];
        fprs.assign(number_expression_thresholds, fpr);
    }
    else if (fprs.size() != number_expression_thresholds)
    {
        throw std::invalid_argument{"Length of false positive rates for IBFs is not equal to length of expression "
                                    "thresholds."};
    }
}

// Estimate the file size for every expression level, necessary when samplewise=false, because then it is completly
// unclear how many minimisers are to store per file.
void get_filsize_per_expression_level(std::filesystem::path const & filename,
                                      uint8_t const number_expression_thresholds,
                                      std::vector<uint16_t> const & expression_thresholds,
                                      std::vector<uint64_t> & sizes,
                                      robin_hood::unordered_set<uint64_t> const & genome,
                                      bool all = true)
{
    std::ifstream fin{filename, std::ios::binary};

    // Skip the first 22 bytes:
    //   8 num_of_minimisers
    //   1 cutoff
    //   1 args.k
    //   4 args.w_size
    //   8 args.s
    fin.ignore(22);

    bool ungapped;
    read_stream(fin, ungapped);
    if (!ungapped)
    {
        fin.ignore(8); // args.shape
    }

    uint64_t minimiser;
    uint16_t minimiser_count;
    sizes.assign(number_expression_thresholds, 0);
    auto const expression_begin_it = expression_thresholds.begin();

    while (read_stream(fin, minimiser))
    {
        read_stream(fin, minimiser_count);
        if (all || genome.contains(minimiser))
        {
            // Find the level with the smallest greater value than the minimiser occurrence, in the level before that the
            // minimiser is going to be stored.
            auto p = std::ranges::upper_bound(expression_thresholds, minimiser_count);
            if (p != expression_begin_it)
                ++sizes[std::ranges::distance(expression_begin_it, p) - 1];
        }
    }
}

// Actual ibf construction
template <bool samplewise, bool minimiser_files_given = true>
void ibf_helper(std::vector<std::filesystem::path> const & minimiser_files,
                std::vector<double> const & fprs,
                estimate_ibf_arguments & ibf_args,
                std::vector<uint8_t> & cutoffs,
                size_t num_hash = 1,
                std::filesystem::path const & expression_by_genome_file = "",
                minimiser_file_input_arguments const & minimiser_args = {})
{
    size_t const num_files = [&]() constexpr
    {
        if constexpr (minimiser_files_given)
            return minimiser_files.size();
        else
            return minimiser_args.samples.size();
    }();

    std::vector<std::vector<uint16_t>> expressions = [&]()
    {
        std::vector<std::vector<uint16_t>> result;
        if constexpr (samplewise)
            result.resize(num_files, std::vector<uint16_t>(ibf_args.number_expression_thresholds));
        return result;
    }();

    std::vector<std::vector<uint64_t>> sizes(num_files);
    std::vector<uint64_t> sizes_ibf{};
    std::vector<std::vector<uint64_t>> counts_per_level(num_files,
                                                        std::vector<uint64_t>(ibf_args.number_expression_thresholds));

    bool const calculate_cutoffs = cutoffs.empty();

    robin_hood::unordered_set<uint64_t> include_set_table; // Storage for minimisers in include file
    robin_hood::unordered_set<uint64_t> exclude_set_table; // Storage for minimisers in exclude file

    if constexpr (!minimiser_files_given)
    {
        if (minimiser_args.include_file != "")
            get_include_set_table(ibf_args, minimiser_args.include_file, include_set_table);
        if (minimiser_args.exclude_file != "")
            get_include_set_table(ibf_args, minimiser_args.exclude_file, exclude_set_table);
    }

    // I/O inside OpenMP region
    if (minimiser_args.ram_friendly)
    {
        seqan3::contrib::bgzf_thread_count = ibf_args.threads;
        omp_set_num_threads(1);
    }
    else
    {
        seqan3::contrib::bgzf_thread_count = 1u;
        omp_set_num_threads(ibf_args.threads);
    }

    // size_t const chunk_size = std::clamp<size_t>(std::bit_ceil(num_files / ibf_args.threads), 8u, 64u);

    // If expression_thresholds should only be depending on minimsers in a certain genome file, genome is created.
    robin_hood::unordered_set<uint64_t> genome{};
    if (expression_by_genome_file != "")
        get_include_set_table(ibf_args, expression_by_genome_file, genome);
    bool const expression_by_genome = (expression_by_genome_file == "");

    // Get expression levels and sizes
    for (size_t i = 0; i < num_files; ++i)
    {
        // Content depends on `minimiser_files_given`:
        // * `false`: filesize
        // * `true` : number of minimisers
        uint64_t filesize{};

        if constexpr (minimiser_files_given)
        {
            uint8_t cutoff;
            read_binary_start(ibf_args, minimiser_files[i], filesize, cutoff);
            cutoffs.push_back(cutoff);
        }
        else
        {
            // Estimate sizes on filesize, assuming every byte translates to one letter (which is obiously not true,
            // because ids contain letters as well), so size might be overestimated. TODO: Find a better estimation!
            size_t const file_iterator =
                std::reduce(minimiser_args.samples.begin(), minimiser_args.samples.begin() + i);

            // Determine cutoffs
            if (calculate_cutoffs)
                cutoffs.push_back(calculate_cutoff(minimiser_files[file_iterator], minimiser_args.samples[i]));

            bool const is_compressed = minimiser_files[file_iterator].extension() == ".gz"
                                    || minimiser_files[file_iterator].extension() == ".bgzf"
                                    || minimiser_files[file_iterator].extension() == ".bz2";
            bool const is_fasta = is_compressed ? check_for_fasta_format(seqan3::format_fasta::file_extensions,
                                                                         minimiser_files[file_iterator].stem())
                                                : check_for_fasta_format(seqan3::format_fasta::file_extensions,
                                                                         minimiser_files[file_iterator].extension());
            filesize = std::filesystem::file_size(minimiser_files[file_iterator]) * minimiser_args.samples[i]
                     * (is_fasta ? 2 : 1) / (is_compressed ? 1 : 3);
            filesize /= ((cutoffs[i] + 1) * (is_fasta ? 1 : 2));
            // ^^ why divide? --> estimate the number of minimisers
        }

        // If set_expression_thresholds_samplewise is not set the expressions as determined by the first file are used for
        // all files.
        if constexpr (samplewise)
        {
            uint64_t diff{1};
            for (int c = 0; c < ibf_args.number_expression_thresholds - 1; c++)
            {
                diff = diff * 2;
                sizes[i].push_back(seqan::hibf::divide_and_ceil(filesize, diff));
            }
            sizes[i].push_back(seqan::hibf::divide_and_ceil(filesize, diff));
        }
        else if constexpr (minimiser_files_given)
        {
            get_filsize_per_expression_level(minimiser_files[i],
                                             ibf_args.number_expression_thresholds,
                                             ibf_args.expression_thresholds,
                                             sizes[i],
                                             genome,
                                             expression_by_genome);
        }
        else
        {
            float diff{1};
            for (int c = 0; c < ibf_args.number_expression_thresholds - 1; c++)
            {
                diff = ibf_args.expression_thresholds[c + 1] / static_cast<float>(ibf_args.expression_thresholds[c]);
                sizes[i].push_back(std::ceil(filesize / diff));
            }
            sizes[i].push_back(std::ceil(filesize / diff));
        }
    }

    size_t const total_bins = num_files * ibf_args.number_expression_thresholds;

    uint64_t max_elements = 0;
    uint64_t total_elements = 0;

    for (size_t i = 0; i < num_files; ++i)
    {
        for (size_t j = 0; j < ibf_args.number_expression_thresholds; ++j)
        {
            max_elements = std::max(max_elements, sizes[i][j]);
            total_elements += sizes[i][j];
        }
    }

    // Check if we have valid data
    if (max_elements == 0)
    {
        throw std::invalid_argument{"No minimizers found for any threshold level."};
    }

    // Use average elements per bin for sizing
    uint64_t avg_elements = std::max(1UL, total_elements / total_bins);
    uint64_t elements_for_sizing = avg_elements;
    // Alternatively: Use maximum elements for sizing:
    // uint64_t elements_for_sizing = std::max(max_elements, avg_elements);

    // Calculate IBF bin size using the first FPR
    double const target_fpr = fprs[0];

    // m = -hn/ln(1-p^(1/h))
    double const log_p_over_h = std::log(target_fpr) / static_cast<double>(num_hash);
    double const one_minus_p_pow_inv_h = 1.0 - std::exp(log_p_over_h);

    if (one_minus_p_pow_inv_h <= 0.0)
    {
        throw std::invalid_argument{"Invalid FPR value - results in invalid bloom filter parameters."};
    }

    double const denominator = std::log(one_minus_p_pow_inv_h);
    double const numerator = -static_cast<double>(elements_for_sizing * num_hash);
    uint64_t const ibf_bin_size = static_cast<uint64_t>(std::ceil(numerator / denominator));

    // Create single IBF
    seqan::hibf::interleaved_bloom_filter single_ibf(seqan::hibf::bin_count{total_bins},
                                                     seqan::hibf::bin_size{ibf_bin_size},
                                                     seqan::hibf::hash_function_count{num_hash});

    // Initialize sizes_ibf for compatibility with FPR calculation
    sizes_ibf.assign(ibf_args.number_expression_thresholds, ibf_bin_size);

    // Collect minimizer insertions from all threads
    std::vector<std::vector<std::pair<uint64_t, std::vector<size_t>>>> file_insertions(num_files);

    // Sequentially process each file TODO: implement parallel processing
    for (size_t i = 0; i < num_files; i++)
    {
        robin_hood::unordered_node_map<uint64_t, uint16_t> hash_table{};
        robin_hood::unordered_node_map<uint64_t, uint8_t> cutoff_table;
        std::vector<uint16_t> expression_thresholds;

        // Fill hash table with minimisers.
        if constexpr (minimiser_files_given)
        {
            read_binary(minimiser_files[i], hash_table);
        }
        else
        {
            size_t const file_iterator =
                std::reduce(minimiser_args.samples.begin(), minimiser_args.samples.begin() + i);
            for (size_t f = 0; f < minimiser_args.samples[i]; f++)
            {
                sequence_file_t fin{minimiser_files[file_iterator + f]};
                if (minimiser_args.ram_friendly)
                    fill_hash_table_parallel(ibf_args,
                                             fin,
                                             hash_table,
                                             cutoff_table,
                                             include_set_table,
                                             exclude_set_table,
                                             (minimiser_args.include_file != ""),
                                             cutoffs[i]);
                else
                    fill_hash_table(ibf_args,
                                    fin,
                                    hash_table,
                                    cutoff_table,
                                    include_set_table,
                                    exclude_set_table,
                                    (minimiser_args.include_file != ""),
                                    cutoffs[i]);
            }
            cutoff_table.clear();
        }

        // Handle samplewise expression thresholds
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

        // Collect insertions for this file
        std::vector<std::pair<uint64_t, std::vector<size_t>>> local_insertions;

        for (auto && [hash, occurrence] : hash_table)
        {
            std::vector<size_t> target_bins;
            for (size_t const j : std::views::iota(0u, ibf_args.number_expression_thresholds) | std::views::reverse)
            {
                uint16_t const threshold = [&]()
                {
                    if constexpr (samplewise)
                        return expressions[i][j];
                    else
                        return ibf_args.expression_thresholds[j];
                }();

                if (occurrence >= threshold)
                {
                    size_t bin_index = j * num_files + i;
                    target_bins.push_back(bin_index);
                    counts_per_level[i][j]++;
                    break;
                }
            }
            if (!target_bins.empty())
            {
                local_insertions.emplace_back(hash, target_bins);
            }
        }
        file_insertions[i] = std::move(local_insertions);
    }

    // Insert all collected minimizers into the IBF (single-threaded)
    for (auto const & file_data : file_insertions)
    {
        for (auto const & [hash, bins] : file_data)
        {
            for (size_t bin_idx : bins)
            {
                single_ibf.emplace(hash, seqan::hibf::bin_index{bin_idx});
            }
        }
    }

    // Store the single IBF
    std::filesystem::path const filename = filenames::ibf(ibf_args.path_out, samplewise, 0, ibf_args);
    store_ibf(single_ibf, filename);

    // Store all expression thresholds per level.
    if constexpr (samplewise)
    {
        std::ofstream outfile{filenames::levels(ibf_args.path_out)};
        for (unsigned j = 0; j < ibf_args.number_expression_thresholds; j++)
        {
            for (size_t i = 0; i < num_files; i++)
                outfile << expressions[i][j] << " ";
            outfile << "\n";
        }
        outfile << "/\n";
    }

    // Store FPR information
    std::ofstream outfile{filenames::fprs(ibf_args.path_out)};
    for (unsigned j = 0; j < ibf_args.number_expression_thresholds; j++)
    {
        for (size_t i = 0; i < num_files; i++)
        {
            // Calculate actual FPR based on IBF parameters
            if (counts_per_level[i][j] > 0 && sizes_ibf[j] > 0)
            {
                double const exp_arg = (num_hash * counts_per_level[i][j]) / static_cast<double>(sizes_ibf[j]);
                double const log_arg = 1.0 - std::exp(-exp_arg);
                double const fpr = (log_arg > 0) ? std::exp(num_hash * std::log(log_arg)) : 1.0;
                outfile << fpr << " ";
            }
            else
            {
                outfile << "0.0 ";
            }
        }
        outfile << "\n";
    }
    outfile << "/\n";
}

// Create ibfs
std::vector<uint16_t> ibf(std::vector<std::filesystem::path> const & sequence_files,
                          estimate_ibf_arguments & ibf_args,
                          minimiser_file_input_arguments & minimiser_args,
                          std::vector<double> & fpr,
                          std::vector<uint8_t> & cutoffs,
                          std::filesystem::path const & expression_by_genome_file,
                          size_t num_hash)
{
    // Declarations
    robin_hood::unordered_node_map<uint64_t, uint16_t> hash_table{}; // Storage for minimisers

    check_cutoffs_samples(sequence_files, minimiser_args.paired, minimiser_args.samples, cutoffs);

    check_expression(ibf_args.expression_thresholds, ibf_args.number_expression_thresholds, expression_by_genome_file);
    check_fpr(ibf_args.number_expression_thresholds, fpr);

    ibf_args.samplewise = (ibf_args.expression_thresholds.size() == 0);

    // Store experiment names
    if (minimiser_args.experiment_names)
    {
        std::ofstream outfile{filenames::stored(ibf_args.path_out)};
        for (size_t i = 0; i < minimiser_args.samples.size(); i++)
        {
            outfile << sequence_files[std::reduce(minimiser_args.samples.begin(), minimiser_args.samples.begin() + i)]
                    << "\n";
        }
    }

    if (ibf_args.samplewise)
        ibf_helper<true, false>(sequence_files,
                                fpr,
                                ibf_args,
                                cutoffs,
                                num_hash,
                                expression_by_genome_file,
                                minimiser_args);
    else
        ibf_helper<false, false>(sequence_files,
                                 fpr,
                                 ibf_args,
                                 cutoffs,
                                 num_hash,
                                 expression_by_genome_file,
                                 minimiser_args);

    store_args(ibf_args, filenames::data(ibf_args.path_out));

    return ibf_args.expression_thresholds;
}

// Create ibfs based on the minimiser file
std::vector<uint16_t> ibf(std::vector<std::filesystem::path> const & minimiser_files,
                          estimate_ibf_arguments & ibf_args,
                          std::vector<double> & fpr,
                          std::filesystem::path const & expression_by_genome_file,
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

    store_args(ibf_args, filenames::data(ibf_args.path_out));

    return ibf_args.expression_thresholds;
}
