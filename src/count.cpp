// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "count.hpp"

#include "misc/fill_hash_table.hpp"
#include "misc/stream.hpp"

void count_genome(minimiser_arguments const & args,
                  std::filesystem::path include_file,
                  std::filesystem::path exclude_file)
{
    robin_hood::unordered_set<uint64_t> include_set_table{};
    robin_hood::unordered_set<uint64_t> exclude_set_table{};

    if (exclude_file != "")
    {
        sequence_file_t fin{exclude_file};
        for (auto & [seq] : fin)
        {
            if (seq.size() >= args.w_size.get())
            {
                for (auto && minHash : seqan3::views::minimiser_hash(seq, args.shape, args.w_size, args.s))
                    exclude_set_table.insert(minHash);
            }
        }
    }

    sequence_file_t fin{include_file};
    for (auto & [seq] : fin)
    {
        if (seq.size() >= args.w_size.get())
        {
            for (auto && minHash : seqan3::views::minimiser_hash(seq, args.shape, args.w_size, args.s))
            {
                if (!(exclude_set_table.contains(minHash)))
                    include_set_table.insert(minHash);
            }
        }
    }

    // Write minimiser to binary
    std::ofstream outfile{std::string{args.path_out} + std::string{include_file.stem()} + ".genome", std::ios::binary};

    for (auto && hash : include_set_table)
    {
        outfile.write(reinterpret_cast<char const *>(&hash), sizeof(hash));
    }
}

void count(minimiser_arguments const & args,
           std::vector<std::filesystem::path> sequence_files,
           std::filesystem::path include_file,
           std::filesystem::path genome_file,
           bool paired)
{
    robin_hood::unordered_node_map<uint64_t, uint16_t> hash_table{};
    // Create a smaller cutoff table to save RAM, this cutoff table is only used for constructing the hash table
    // and afterwards discarded.
    robin_hood::unordered_node_map<uint64_t, uint8_t> cutoff_table;
    robin_hood::unordered_set<uint64_t> include_set_table{};
    robin_hood::unordered_set<uint64_t> exclude_set_table{};
    std::vector<uint64_t> counter{};
    uint64_t expression{};

    // Read minimiser from binary
    {
        std::ifstream fin{genome_file, std::ios::binary};

        uint64_t minimiser;
        while (read_stream(fin, minimiser))
        {
            include_set_table.insert(minimiser);
        }
    }

    for (size_t i = 0; i < sequence_files.size(); ++i)
    {
        if (paired)
        {
            sequence_file_t fin{sequence_files[i]};
            fill_hash_table(args, fin, hash_table, cutoff_table, include_set_table, exclude_set_table, true);
            ++i;
            fin = sequence_files[i];
            fill_hash_table(args, fin, hash_table, cutoff_table, include_set_table, exclude_set_table, true);
        }
        else
        {
            sequence_file_t fin{sequence_files[i]};
            fill_hash_table(args, fin, hash_table, cutoff_table, include_set_table, exclude_set_table, true);
        }
        cutoff_table.clear();

        std::ofstream outfile{std::string{args.path_out} + std::string{sequence_files[i].stem()} + ".count.out"};
        sequence_file_with_id_t fin{include_file};

        for (auto & [id, seq] : fin)
        {
            if (seq.size() >= args.w_size.get())
            {
                for (auto && minHash : seqan3::views::minimiser_hash(seq, args.shape, args.w_size, args.s))
                    counter.push_back(hash_table[minHash]);
                std::ranges::nth_element(counter, counter.begin() + counter.size() / 2);
                expression = counter[counter.size() / 2];
                outfile << id << "\t" << expression << "\n";
                counter.clear();
            }
        }
        hash_table.clear();
    }
}
