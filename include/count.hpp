// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <filesystem>
#include <vector>

#include "shared.hpp"

/*!\brief Creates a set of minimizers to ignore, which should be used as an input to count.
* \param args               The minimiser arguments to use (seed, shape, window size).
* \param include_file        A file containing the transcripts which expression values should be determined.
* \param exclude_file       A file containing minimizers which should be ignored.
*/
void count_genome(minimiser_arguments const & args,
                  std::filesystem::path const & include_file,
                  std::filesystem::path const & exclude_file);

/*!\brief Get the concrete expression values (= median of all counts of one transcript) for given experiments.
*         This function can be used to estimate how good the median approach can be, if all count values are available.
* \param args               The minimiser arguments to use (seed, shape, window size).
* \param sequence_files     The sequence files, which contains the reads.
* \param include_file       A file containing the transcripts which expression values should be determined.
* \param genome_file        A "*.genome" file constructed with the command genome.
* \param paired             Flag to indicate if input data is paired or not.
*/
void count(minimiser_arguments const & args,
           std::vector<std::filesystem::path> const & sequence_files,
           std::filesystem::path const & include_file,
           std::filesystem::path const & genome_file,
           bool paired);
