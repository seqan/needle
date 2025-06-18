// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "insert.hpp"

#include "misc/calculate_cutoff.hpp"
#include "misc/check_cutoffs_samples.hpp"
#include "misc/check_for_fasta_format.hpp"
#include "misc/filenames.hpp"
#include "misc/fill_hash_table.hpp"
#include "misc/get_expression_thresholds.hpp"
#include "misc/get_include_set_table.hpp"
#include "misc/read_levels.hpp"
#include "misc/stream.hpp"

// Actual insertion
template <bool samplewise, bool minimiser_files_given = true>
void insert_helper(std::vector<std::filesystem::path> const & minimiser_files,
                   estimate_ibf_arguments & ibf_args,
                   std::filesystem::path path_in,
                   std::vector<uint8_t> & cutoffs,
                   std::filesystem::path const & expression_by_genome_file = "",
                   minimiser_file_input_arguments const & minimiser_args = {})
{
    size_t old_bin_number{};
    size_t new_bin_number{};
    size_t num_hash_functions{};

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

    std::vector<uint64_t> sizes{};
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

    size_t const chunk_size = std::clamp<size_t>(std::bit_ceil(num_files / ibf_args.threads), 8u, 64u);

    // If expression_thresholds should only be depending on minimsers in a certain genome file, genome is created.
    robin_hood::unordered_set<uint64_t> genome{};
    if (expression_by_genome_file != "")
        get_include_set_table(ibf_args, expression_by_genome_file, genome);
    bool const expression_by_genome = (expression_by_genome_file == "");

    // Check, if there are deleted bins
    std::vector<uint64_t> pos_insert{};
    if (std::filesystem::path const deleted_files_path = filenames::deleted(path_in);
        std::filesystem::exists(deleted_files_path))
    {
        std::ifstream fin{deleted_files_path};
        uint64_t number;

        while (fin >> number)
        {
            pos_insert.push_back(number);
        }
    }
    size_t num_deleted = pos_insert.size();

    // Adjust IBFs
    std::vector<seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>> ibfs;
    for (size_t j = 0; j < ibf_args.number_expression_thresholds; j++)
    {
        seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf_load{};
        load_ibf(ibf_load, filenames::ibf(path_in, samplewise, j, ibf_args));

        old_bin_number = ibf_load.bin_count();
        new_bin_number = old_bin_number + num_files - num_deleted;
        if (new_bin_number > old_bin_number)
            ibf_load.increase_bin_number_to(seqan3::bin_count{new_bin_number});

        num_hash_functions = ibf_load.hash_function_count();
        sizes.push_back(ibf_load.bin_size());
        ibfs.push_back(std::move(ibf_load));
    }

    for (size_t j = old_bin_number; j < new_bin_number; j++)
        pos_insert.push_back(j);

// Add minimisers to ibf
#pragma omp parallel for schedule(dynamic, chunk_size)
    for (size_t i = 0; i < num_files; i++)
    {
        robin_hood::unordered_node_map<uint64_t, uint16_t> hash_table{}; // Storage for minimisers
        // Create a smaller cutoff table to save RAM, this cutoff table is only used for constructing the hash table
        // and afterwards discarded.
        robin_hood::unordered_node_map<uint64_t, uint8_t> cutoff_table;
        std::vector<uint16_t> expression_thresholds;

        // Fill hash table with minimisers.
        if constexpr (minimiser_files_given)
        {
            read_binary(minimiser_files[i], hash_table);
            uint8_t cutoff;
            uint64_t tmp_var{};
            read_binary_start(ibf_args, minimiser_files[i], tmp_var, cutoff);
            cutoffs.push_back(cutoff);
        }
        else
        {
            // Estimate sizes on filesize, assuming every byte translates to one letter (which is obiously not true,
            // because ids contain letters as well), so size might be overestimated. TODO: Find a better estimation!
            size_t const file_iterator =
                std::reduce(minimiser_args.samples.begin(), minimiser_args.samples.begin() + i);
            uint64_t filesize{};
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
            filesize = filesize / ((cutoffs[i] + 1) * (is_fasta ? 1 : 2));

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

        // If set_expression_thresholds_samplewise is not set the expressions as determined by the first file are used for
        // all files.
        if constexpr (samplewise)
        {
            std::vector<uint64_t> sizes_tmp{};
            get_expression_thresholds(ibf_args.number_expression_thresholds,
                                      hash_table,
                                      expression_thresholds,
                                      sizes_tmp,
                                      genome,
                                      cutoffs[i],
                                      expression_by_genome);
            expressions[i] = expression_thresholds;
            sizes_tmp.clear();
        }

        // Every minimiser is stored in IBF, if it occurence is greater than or equal to the expression level
        for (auto && elem : hash_table)
        {
            for (size_t const j : std::views::iota(0u, ibf_args.number_expression_thresholds) | std::views::reverse)
            {
                uint16_t const threshold = [&]()
                {
                    if constexpr (samplewise)
                        return expressions[i][j];
                    else
                        return ibf_args.expression_thresholds[j];
                }();

                if (elem.second >= threshold)
                {
                    ibfs[j].emplace(elem.first, seqan3::bin_index{pos_insert[i]});
                    counts_per_level[i][j]++;
                    break;
                }
            }
        }
    }

    // Store IBFs
    for (unsigned i = 0; i < ibf_args.number_expression_thresholds; i++)
    {
        std::filesystem::path const filename = filenames::ibf(ibf_args.path_out, samplewise, i, ibf_args);

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
    if constexpr (samplewise)
    {
        std::vector<std::vector<uint16_t>> expressions_prev{};
        read_levels<uint16_t>(expressions_prev, filenames::levels(path_in));

        std::ofstream outfile{filenames::levels(ibf_args.path_out)};
        for (unsigned j = 0; j < ibf_args.number_expression_thresholds; j++)
        {
            int exp_i = 0;
            for (unsigned i = 0; i < new_bin_number; i++)
            {
                if (std::ranges::find(pos_insert, i) != pos_insert.end())
                {
                    outfile << expressions[exp_i][j] << " ";
                    exp_i++;
                }
                else
                {
                    outfile << expressions_prev[j][i] << " ";
                }
            }
            outfile << "\n";
        }
        outfile << "/\n";
    }

    // Store all fprs per level.
    std::vector<std::vector<double>> fprs_prev{};
    read_levels<double>(fprs_prev, filenames::fprs(path_in));

    std::ofstream outfile{filenames::fprs(ibf_args.path_out)};
    for (unsigned j = 0; j < ibf_args.number_expression_thresholds; j++)
    {
        size_t exp_i = 0;
        for (size_t i = 0; i < new_bin_number; i++)
        {
            if (std::ranges::find(pos_insert, i) != pos_insert.end())
            {
                // m = -hn/ln(1-p^(1/h))
                double const exp_arg =
                    (num_hash_functions * counts_per_level[exp_i][j]) / static_cast<double>(sizes[j]);
                double const log_arg = 1.0 - std::exp(-exp_arg);
                double const fpr = std::exp(num_hash_functions * std::log(log_arg));
                outfile << fpr << " ";
                exp_i++;
            }
            else
            {
                outfile << fprs_prev[j][i] << " ";
            }
        }
        outfile << "\n";
    }
    outfile << "/\n";
}

// Insert into ibfs
std::vector<uint16_t> insert(std::vector<std::filesystem::path> const & sequence_files,
                             estimate_ibf_arguments & ibf_args,
                             minimiser_file_input_arguments & minimiser_args,
                             std::vector<uint8_t> & cutoffs,
                             std::filesystem::path const & expression_by_genome_file,
                             std::filesystem::path const & path_in,
                             bool samplewise)
{
    // Declarations
    robin_hood::unordered_node_map<uint64_t, uint16_t> hash_table{}; // Storage for minimisers

    check_cutoffs_samples(sequence_files, minimiser_args.paired, minimiser_args.samples, cutoffs);
    load_args(ibf_args, filenames::data(path_in));

    if (samplewise)
        insert_helper<true, false>(sequence_files,
                                   ibf_args,
                                   path_in,
                                   cutoffs,
                                   expression_by_genome_file,
                                   minimiser_args);
    else
        insert_helper<false, false>(sequence_files,
                                    ibf_args,
                                    path_in,
                                    cutoffs,
                                    expression_by_genome_file,
                                    minimiser_args);

    store_args(ibf_args, filenames::data(ibf_args.path_out));
    return ibf_args.expression_thresholds;
}

// Insert into ibfs based on the minimiser file
std::vector<uint16_t> insert(std::vector<std::filesystem::path> const & minimiser_files,
                             estimate_ibf_arguments & ibf_args,
                             std::filesystem::path const & expression_by_genome_file,
                             std::filesystem::path const & path_in,
                             bool samplewise)
{
    std::vector<uint8_t> cutoffs{};
    load_args(ibf_args, filenames::data(path_in));
    if (samplewise)
        insert_helper<true>(minimiser_files, ibf_args, path_in, cutoffs, expression_by_genome_file);
    else
        insert_helper<false>(minimiser_files, ibf_args, path_in, cutoffs, expression_by_genome_file);

    store_args(ibf_args, filenames::data(ibf_args.path_out));

    return ibf_args.expression_thresholds;
}
