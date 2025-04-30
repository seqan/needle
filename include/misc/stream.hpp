// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstdint>
#include <filesystem>
#include <iostream>
#include <robin_hood.h>

#include "shared.hpp"

template <typename t>
inline auto & read_stream(std::ifstream & stream, t & val)
{
    return stream.read(reinterpret_cast<char *>(std::addressof(val)), sizeof(val));
}

template <typename t>
inline auto & write_stream(std::ofstream & stream, t const & val)
{
    return stream.write(reinterpret_cast<char const *>(std::addressof(val)), sizeof(val));
}

/*!\brief Reads a binary file that needle minimiser creates.
* \param filename           The filename of the binary file.
* \param hash_table         The hash table to store minimisers into.

*/
void read_binary(std::filesystem::path filename, robin_hood::unordered_node_map<uint64_t, uint16_t> & hash_table);

/*!\brief Reads the beginning of a binary file that needle minimiser creates.
* \param args               Min arguments.
* \param filename           The filename of the binary file.
* \param num_of_minimisers  Variable, where to number of minimisers should be stored.
* \param cutoff             cutoff value.
*/
void read_binary_start(minimiser_arguments & args,
                       std::filesystem::path filename,
                       uint64_t & num_of_minimisers,
                       uint8_t & cutoff);
