// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "minimiser.hpp"

#include <fstream>
#include <numeric>
#include <omp.h>
#include <utility>

#include "misc/calculate_cutoff.hpp"
#include "misc/check_cutoffs_samples.hpp"
#include "misc/fill_hash_table.hpp"
#include "misc/get_expression_thresholds.hpp"
#include "misc/get_include_set_table.hpp"
#include "misc/stream.hpp"

// Add result struct
struct MinResult
{
    uint64_t count{};
    uint8_t cutoff{};
    std::vector<uint16_t> expression_thresholds;
    std::vector<uint64_t> sizes;
};

// Actuall minimiser calculation
template <bool parallel = false>
MinResult calculate_minimiser(std::vector<std::filesystem::path> const & sequence_files,
                              robin_hood::unordered_set<uint64_t> const & include_set_table,
                              robin_hood::unordered_set<uint64_t> const & exclude_set_table,
                              minimiser_arguments const & args,
                              minimiser_file_input_arguments const & minimiser_args,
                              unsigned const i,
                              std::vector<uint8_t> & cutoffs,
                              uint8_t const number_expression_thresholds,
                              bool const write_counts)
{
    robin_hood::unordered_node_map<uint64_t, uint16_t> hash_table{}; // Storage for minimisers
    uint8_t cutoff{0};

    // Create a smaller cutoff table to save RAM, this cutoff table is only used for constructing the hash table
    // and afterwards discarded.
    robin_hood::unordered_node_map<uint64_t, uint8_t> cutoff_table;
    size_t const file_iterator = std::reduce(minimiser_args.samples.begin(), minimiser_args.samples.begin() + i);

    bool const calculate_cutoffs = cutoffs.empty();

    if (calculate_cutoffs)
        cutoff = calculate_cutoff(sequence_files[file_iterator], minimiser_args.samples[i]);
    else
        cutoff = cutoffs[i];

    // Fill hash_table with minimisers.
    for (size_t f = 0; f < minimiser_args.samples[i]; f++)
    {
        sequence_file_t fin{sequence_files[file_iterator + f]};
        if constexpr (parallel)
        {
            fill_hash_table_parallel(args,
                                     fin,
                                     hash_table,
                                     cutoff_table,
                                     include_set_table,
                                     exclude_set_table,
                                     (minimiser_args.include_file != ""),
                                     cutoff);
        }
        else
        {
            fill_hash_table(args,
                            fin,
                            hash_table,
                            cutoff_table,
                            include_set_table,
                            exclude_set_table,
                            (minimiser_args.include_file != ""),
                            cutoff);
        }
    }
    cutoff_table.clear();

    // Write minimiser and their counts to binary
    std::ofstream outfile{std::string{args.path_out} + std::string{sequence_files[file_iterator].stem()} + ".minimiser",
                          std::ios::binary};
    auto hash_size = hash_table.size();

    write_stream(outfile, hash_size);
    write_stream(outfile, cutoff);
    write_stream(outfile, args.k);
    write_stream(outfile, args.w_size);
    write_stream(outfile, args.s);
    bool ungapped = args.shape.all();
    write_stream(outfile, ungapped);

    if (!ungapped)
        write_stream(outfile, args.shape);

    for (auto && hash : hash_table)
    {
        write_stream(outfile, hash.first);
        write_stream(outfile, hash.second);
    }

    uint64_t count = hash_table.size();

    MinResult res;
    res.count = count;
    res.cutoff = cutoff;

    // If insert counts requested, compute expression thresholds (sizes are computed).
    if (write_counts)
        get_expression_thresholds(number_expression_thresholds,
                                  hash_table,
                                  res.expression_thresholds,
                                  res.sizes,
                                  robin_hood::unordered_set<uint64_t>{},
                                  cutoff,
                                  true);

    return res;
}

void minimiser(std::vector<std::filesystem::path> const & sequence_files,
               minimiser_arguments const & args,
               minimiser_file_input_arguments & minimiser_args,
               std::vector<uint8_t> & cutoffs)
{
    // Declarations
    robin_hood::unordered_set<uint64_t> include_set_table{}; // Storage for minimisers in include file
    robin_hood::unordered_set<uint64_t> exclude_set_table{}; // Storage for minimisers in exclude file

    check_cutoffs_samples(sequence_files, minimiser_args.paired, minimiser_args.samples, cutoffs);

    if (minimiser_args.include_file != "")
        get_include_set_table(args, minimiser_args.include_file, include_set_table);
    if (minimiser_args.exclude_file != "")
        get_include_set_table(args, minimiser_args.exclude_file, exclude_set_table);

    size_t const chunk_size = std::clamp<size_t>(std::bit_ceil(minimiser_args.samples.size() / args.threads), 1u, 64u);

    // Prepare container for per-user-bin counts
    std::vector<uint64_t> insert_counts(minimiser_args.samples.size(), 0);
    std::vector<std::vector<uint16_t>> all_thresholds;
    std::vector<std::vector<uint64_t>> all_sizes;
    if (minimiser_args.write_counts)
    {
        all_thresholds.resize(minimiser_args.samples.size());
        all_sizes.resize(minimiser_args.samples.size());
    }

    // Add minimisers to ibf
    if (minimiser_args.ram_friendly)
    {
        seqan3::contrib::bgzf_thread_count = args.threads;
        for (size_t i = 0; i < minimiser_args.samples.size(); i++)
        {
            MinResult r = calculate_minimiser<true>(sequence_files,
                                                    include_set_table,
                                                    exclude_set_table,
                                                    args,
                                                    minimiser_args,
                                                    i,
                                                    cutoffs,
                                                    minimiser_args.number_expression_thresholds,
                                                    minimiser_args.write_counts);
            insert_counts[i] = r.count;
            if (cutoffs.size() <= i)
                cutoffs.resize(minimiser_args.samples.size());
            cutoffs[i] = r.cutoff;
            if (minimiser_args.write_counts)
            {
                all_thresholds[i] = std::move(r.expression_thresholds);
                all_sizes[i] = std::move(r.sizes);
            }
        }
    }
    else
    {
        // I/O inside OpenMP region
        seqan3::contrib::bgzf_thread_count = 1u;
        omp_set_num_threads(args.threads);

        std::vector<MinResult> results(minimiser_args.samples.size());

#pragma omp parallel for schedule(dynamic, chunk_size)
        for (size_t i = 0; i < minimiser_args.samples.size(); i++)
        {
            results[i] = calculate_minimiser<false>(sequence_files,
                                                    include_set_table,
                                                    exclude_set_table,
                                                    args,
                                                    minimiser_args,
                                                    i,
                                                    cutoffs,
                                                    minimiser_args.number_expression_thresholds,
                                                    minimiser_args.write_counts);
        }

        for (size_t i = 0; i < minimiser_args.samples.size(); ++i)
        {
            insert_counts[i] = results[i].count;
            if (cutoffs.size() <= i)
                cutoffs.resize(minimiser_args.samples.size());
            cutoffs[i] = results[i].cutoff;
            if (minimiser_args.write_counts)
            {
                all_thresholds[i] = std::move(results[i].expression_thresholds);
                all_sizes[i] = std::move(results[i].sizes);
            }
        }
    }

    // If --write-counts is set, write thresholds.tsv and per-user-bin count
    if (minimiser_args.write_counts)
    {
        // write thresholds.tsv
        std::filesystem::path thr_path = std::filesystem::path{args.path_out} / "thresholds.tsv";
        std::ofstream thr_out{thr_path};
        thr_out << "file";
        for (uint8_t l = 0; l < minimiser_args.number_expression_thresholds; ++l)
            thr_out << '\t' << "level_" << static_cast<int>(l);
        thr_out << '\n';

        for (size_t i = 0, file_it = 0; i < minimiser_args.samples.size(); ++i)
        {
            thr_out << sequence_files[file_it].filename().string();
            for (uint8_t l = 0; l < minimiser_args.number_expression_thresholds; ++l)
            {
                uint16_t val = 0;
                if (l < all_thresholds[i].size())
                    val = all_thresholds[i][l];
                thr_out << '\t' << val;
            }
            thr_out << '\n';
            file_it += minimiser_args.samples[i];
        }

        // Write insert counts
        // Compute per level number of inserts by threshold positions from get_expression_thresholds
        std::filesystem::path sizes_path = std::filesystem::path{args.path_out} / "thresholds_inserts.tsv";
        std::ofstream sizes_out{sizes_path};

        for (size_t i = 0, file_it = 0; i < minimiser_args.samples.size(); ++i)
        {
            sizes_out << sequence_files[file_it].filename().string();
            uint64_t prev_pos = 0;
            size_t const L = static_cast<size_t>(minimiser_args.number_expression_thresholds);
            for (size_t l = 0; l < L; ++l)
            {
                uint64_t pos = 0;
                if (l < all_sizes[i].size())
                    pos = all_sizes[i][l];

                uint64_t per_level = 0;
                if (l + 1 < L)
                {
                    per_level = (pos > prev_pos) ? (pos - prev_pos) : 0;
                    prev_pos = pos;
                }
                else
                {
                    uint64_t total_inserts = (i < insert_counts.size()) ? insert_counts[i] : 0;
                    per_level = (total_inserts > prev_pos) ? (total_inserts - prev_pos) : 0;
                }

                sizes_out << '\t' << per_level;
            }
            sizes_out << '\n';
            file_it += minimiser_args.samples[i];
        }
    }
}
