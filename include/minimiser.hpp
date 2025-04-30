// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstdint>
#include <filesystem>
#include <vector>

#include "shared.hpp"

/*! \brief Create minimiser and header files.
* \param sequence_files  A vector of sequence file paths.
* \param args            The minimiser arguments to use (seed, shape, window size).
* \param minimiser_args  The minimiser specific arguments to use.
* \param cutoffs         List of cutoffs.
*/
void minimiser(std::vector<std::filesystem::path> const & sequence_files,
               minimiser_arguments const & args,
               minimiser_file_input_arguments & minimiser_args,
               std::vector<uint8_t> & cutoffs);
