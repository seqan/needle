// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "delete.hpp"

// Delete from ibfs
void delete_bin(std::vector<uint64_t> const & delete_files,
                estimate_ibf_arguments & ibf_args,
                std::filesystem::path path_in,
                bool samplewise)
{
    load_args(ibf_args, std::string{path_in} + "IBF_Data");

    std::vector<seqan3::bin_index> bins_to_delete{};
    for (size_t i = 0; i < delete_files.size(); i++)
        bins_to_delete.push_back(seqan3::bin_index{delete_files[i]});

    omp_set_num_threads(ibf_args.threads);

// Delete bins from ibfs
#pragma omp parallel
    for (unsigned i = 0; i < ibf_args.number_expression_thresholds; i++)
    {
        std::filesystem::path filename;
        if (samplewise)
            filename = path_in.string() + "IBF_Level_" + std::to_string(i);
        else
            filename = path_in.string() + "IBF_" + std::to_string(ibf_args.expression_thresholds[i]);

        seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf;
        load_ibf(ibf, filename);
        ibf.clear(bins_to_delete);

        if (samplewise)
            filename = ibf_args.path_out.string() + "IBF_Level_" + std::to_string(i);
        else
            filename = ibf_args.path_out.string() + "IBF_" + std::to_string(ibf_args.expression_thresholds[i]);

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
    std::ofstream outfile{std::string{ibf_args.path_out} + "IBF_Deleted"};
    for (auto & deleted : delete_files)
    {
        outfile << deleted << ",";
    }
    outfile << "\n";
}
