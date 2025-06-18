// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "minimiser.hpp"

#include <numeric>
#include <omp.h>

#include "misc/calculate_cutoff.hpp"
#include "misc/check_cutoffs_samples.hpp"
#include "misc/fill_hash_table.hpp"
#include "misc/get_include_set_table.hpp"
#include "misc/stream.hpp"

// Actuall minimiser calculation
template <bool parallel = false>
void calculate_minimiser(std::vector<std::filesystem::path> const & sequence_files,
                         robin_hood::unordered_set<uint64_t> const & include_set_table,
                         robin_hood::unordered_set<uint64_t> const & exclude_set_table,
                         minimiser_arguments const & args,
                         minimiser_file_input_arguments const & minimiser_args,
                         unsigned const i,
                         std::vector<uint8_t> & cutoffs)
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

    // Add minimisers to ibf
    if (minimiser_args.ram_friendly)
    {
        seqan3::contrib::bgzf_thread_count = args.threads;
        for (size_t i = 0; i < minimiser_args.samples.size(); i++)
        {
            calculate_minimiser<true>(sequence_files,
                                      include_set_table,
                                      exclude_set_table,
                                      args,
                                      minimiser_args,
                                      i,
                                      cutoffs);
        }
    }
    else
    {
        // I/O inside OpenMP region
        seqan3::contrib::bgzf_thread_count = 1u;
        omp_set_num_threads(args.threads);

#pragma omp parallel for schedule(dynamic, chunk_size)
        for (size_t i = 0; i < minimiser_args.samples.size(); i++)
        {
            calculate_minimiser(sequence_files, include_set_table, exclude_set_table, args, minimiser_args, i, cutoffs);
        }
    }
}
