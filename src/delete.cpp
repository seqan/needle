// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "delete.hpp"

#include "misc/filenames.hpp"

// Delete from ibfs
void delete_bin(std::vector<uint64_t> const & delete_files,
                estimate_ibf_arguments & ibf_args,
                std::filesystem::path const & path_in,
                bool samplewise)
{
    load_args(ibf_args, filenames::data(path_in));

    std::vector<seqan3::bin_index> bins_to_delete{};
    bins_to_delete.reserve(delete_files.size());
    for (size_t i = 0; i < delete_files.size(); i++)
        bins_to_delete.push_back(seqan3::bin_index{delete_files[i]});

    omp_set_num_threads(ibf_args.threads);

// Delete bins from ibfs
#pragma omp parallel
    for (unsigned i = 0; i < ibf_args.number_expression_thresholds; i++)
    {
        std::filesystem::path filename = filenames::ibf(path_in, samplewise, i, ibf_args);

        seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf;
        load_ibf(ibf, filename);
        ibf.clear(bins_to_delete);

        filename = filenames::ibf(ibf_args.path_out, samplewise, i, ibf_args);

        if (ibf_args.compressed)
        {
            seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibfc{std::move(ibf)};
            store_ibf(ibfc, filename);
        }
        else
        {
            store_ibf(ibf, filename);
        }
    }

    // Store deleted bins
    std::ofstream outfile{filenames::deleted(ibf_args.path_out)};
    for (auto & deleted : delete_files)
    {
        outfile << deleted << ",";
    }
    outfile << "\n";
}
