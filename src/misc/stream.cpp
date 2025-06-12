// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "misc/stream.hpp"

void read_binary(std::filesystem::path const & filename,
                 robin_hood::unordered_node_map<uint64_t, uint16_t> & hash_table)
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

    while (read_stream(fin, minimiser))
    {
        read_stream(fin, minimiser_count);
        hash_table[minimiser] = minimiser_count;
    }
}

void read_binary_start(minimiser_arguments & args,
                       std::filesystem::path const & filename,
                       uint64_t & num_of_minimisers,
                       uint8_t & cutoff)
{
    std::ifstream fin{filename, std::ios::binary};

    read_stream(fin, num_of_minimisers);
    read_stream(fin, cutoff);
    read_stream(fin, args.k);
    read_stream(fin, args.w_size);
    read_stream(fin, args.s);

    bool ungapped;
    read_stream(fin, ungapped);
    if (ungapped)
        args.shape = seqan3::ungapped{args.k};
    else
        read_stream(fin, args.shape);
}
