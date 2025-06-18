// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstdint>
#include <filesystem>
#include <vector>

#include "shared.hpp"

/*! \brief Delete bins from ibfs
* \param delete_files    A vector of integers specifiying the bins to delete.
* \param ibf_args        The IBF specific arguments to use (bin size, number of hash functions, ...). See
*                        struct ibf_arguments.
* \param path_in         Input directory.
* \param samplewise      True, if expression levels were set beforehand.
*/
void delete_bin(std::vector<uint64_t> const & delete_files,
                estimate_ibf_arguments & ibf_args,
                std::filesystem::path const & path_in,
                bool samplewise);
