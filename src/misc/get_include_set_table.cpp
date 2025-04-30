// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "misc/get_include_set_table.hpp"

// Create set with hashes from the minimisers from an include or exclude file.
void get_include_set_table(minimiser_arguments const & args,
                           std::filesystem::path const include_file,
                           robin_hood::unordered_set<uint64_t> & include_table)
{
    sequence_file_t fin{include_file};
    for (auto & [seq] : fin)
    {
        if (seq.size() >= args.w_size.get())
        {
            for (auto && minHash : seqan3::views::minimiser_hash(seq, args.shape, args.w_size, args.s))
                include_table.insert(minHash);
        }
    }
}
