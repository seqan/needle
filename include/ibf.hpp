// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <filesystem>
#include <iostream>
#include <math.h>
#include <numeric>
#include <string>

#include <seqan3/alphabet/nucleotide/dna4.hpp>

#include "shared.hpp"

/*! \brief Creates IBFs.
 * \param sequence_files  A vector of sequence file paths.
 * \param ibf_args        The IBF specific arguments to use (bin size, number of hash functions, ...). See
 *                        struct ibf_arguments.
 * \param minimiser_args  The minimiser specific arguments to use.
 * \param fpr             The average false positive rate that should be used.
 * \param cutoffs         List of cutoffs.
 * \param expression_by_genome_file File that contains the only minimisers that should be considered for the
 *                                  determination of the expression thresholds.
 * \param num_hash        The number of hash functions to use.
 *  \returns The expression thresholds per experiment.
 */
std::vector<uint16_t> ibf(std::vector<std::filesystem::path> const & sequence_files,
                          estimate_ibf_arguments & ibf_args,
                          minimiser_file_input_arguments & minimiser_args,
                          std::vector<double> & fpr,
                          std::vector<uint8_t> & cutoffs,
                          std::filesystem::path const & expression_by_genome_file = "",
                          size_t num_hash = 1);

/*! \brief Creates IBFs based on the minimiser files
 * \param minimiser_files A vector of minimiser file paths.
 * \param ibf_args        The IBF specific arguments to use (bin size, number of hash functions, ...). See
 *                        struct ibf_arguments.
 * \param fpr             The average false positive rate that should be used.
 * \param expression_by_genome_file File that contains the only minimisers that should be comnsidered for the
 *                                  determination of the expression_thresholds.
 * \param num_hash        The number of hash functions to use.
 *  \returns The expression thresholds per experiment.
 */
std::vector<uint16_t> ibf(std::vector<std::filesystem::path> const & minimiser_files,
                          estimate_ibf_arguments & ibf_args,
                          std::vector<double> & fpr,
                          std::filesystem::path const & expression_by_genome_file = "",
                          size_t num_hash = 1);
