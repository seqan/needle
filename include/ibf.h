// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/needle/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <iostream>
#include <math.h>
#include <numeric>
#include <string>

#include <seqan3/alphabet/container/concatenated_sequences.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <filesystem>

#include "shared.h"

struct minimiser_arguments
{
    std::filesystem::path include_file; // Needs to be defined when only minimisers appearing in this file should be stored
    std::filesystem::path exclude_file; // Needs to be defined when minimisers appearing in this file should NOT be stored
    std::vector<int> samples{}; // Can be used to indicate that sequence files belong to the same experiment
    bool paired = false; // If true, than experiments are seen as paired-end experiments
    bool experiment_names = false; // Flag, if names of experiment should be stored in a txt file
    bool ram_friendly = false;
};

//!\brief Generates a random integer not greater than a given maximum
struct RandomGenerator {
	int maxi;
	RandomGenerator(int max) :
			maxi(max) {
	}

	int operator()() {
		return rand() % maxi;
	}
};

/*!\brief Get the concrete expression values (= median of all counts of one transcript) for given experiments.
*         This function can be used to estimate how good the median approach can be, if all count values are available.
* \param args               The minimiser arguments to use (seed, shape, window size).
* \param sequence_files     The sequence files, which contains the reads.
* \param include_file       A file containing the transcripts which expression values should be determined.
* \param genome_file        A "*.genome" file constructed with the command genome.
* \param paired             Flag to indicate if input data is paired or not.
*/
void count(min_arguments const & args, std::vector<std::filesystem::path> sequence_files, std::filesystem::path include_file,
           std::filesystem::path genome_file, bool paired);

/*!\brief Creates a set of minimizers to ignore, which should be used as an input to count.
* \param args               The minimiser arguments to use (seed, shape, window size).
* \param include_file        A file containing the transcripts which expression values should be determined.
* \param exclude_file       A file containing minimizers which should be ignored.
*/
void count_genome(min_arguments const & args, std::filesystem::path include_file, std::filesystem::path exclude_file);

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
void read_binary_start(min_arguments & args, std::filesystem::path filename, uint64_t & num_of_minimisers, uint8_t & cutoff);

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
std::vector<uint16_t> ibf(std::vector<std::filesystem::path> const & sequence_files, estimate_ibf_arguments & ibf_args,
                          minimiser_arguments & minimiser_args, std::vector<double> & fpr, std::vector<uint8_t> & cutoffs,
                          std::filesystem::path const expression_by_genome_file = "",
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
                          estimate_ibf_arguments & ibf_args, std::vector<double> & fpr,
                          std::filesystem::path const expression_by_genome_file = "",
                          size_t num_hash = 1);

/*! \brief Create minimiser and header files.
* \param sequence_files  A vector of sequence file paths.
* \param args            The minimiser arguments to use (seed, shape, window size).
* \param minimiser_args  The minimiser specific arguments to use.
* \param cutoffs         List of cutoffs.
*/
void minimiser(std::vector<std::filesystem::path> const & sequence_files, min_arguments const & args,
               minimiser_arguments & minimiser_args, std::vector<uint8_t> & cutoffs);

/*! \brief Insert into IBFs.
* \param sequence_files  A vector of sequence file paths.
* \param ibf_args        The IBF specific arguments to use (bin size, number of hash functions, ...). See
*                        struct ibf_arguments.
* \param minimiser_args  The minimiser specific arguments to use.
* \param cutoffs         List of cutoffs.
* \param expression_by_genome_file File that contains the only minimisers that should be considered for the
*                                  determination of the expression thresholds.
* \param path_in         Input directory.
* \param samplewise      True, if expression levels were set beforehand.
*  \returns The expression thresholds per experiment.
*/
std::vector<uint16_t> insert(std::vector<std::filesystem::path> const & sequence_files,
                             estimate_ibf_arguments & ibf_args, minimiser_arguments & minimiser_args,
                             std::vector<uint8_t> & cutoffs,
                             std::filesystem::path const expression_by_genome_file, std::filesystem::path path_in, bool samplewise);

/*! \brief Insert into IBFs based on the minimiser files
* \param minimiser_files A vector of minimiser file paths.
* \param ibf_args        The IBF specific arguments to use (bin size, number of hash functions, ...). See
*                        struct ibf_arguments.
* \param expression_by_genome_file File that contains the only minimisers that should be comnsidered for the
*                                  determination of the expression_thresholds.
* \param path_in         Input directory.
* \param samplewise      True, if expression levels were set beforehand.
*  \returns The expression thresholds per experiment.
*/
std::vector<uint16_t> insert(std::vector<std::filesystem::path> const & minimiser_files,
                             estimate_ibf_arguments & ibf_args,
                             std::filesystem::path const expression_by_genome_file, std::filesystem::path path_in, bool samplewise);

/*! \brief Delete bins from ibfs
* \param delete_files    A vector of integers specifiying the bins to delete.
* \param ibf_args        The IBF specific arguments to use (bin size, number of hash functions, ...). See
*                        struct ibf_arguments.
* \param path_in         Input directory.
* \param samplewise      True, if expression levels were set beforehand.
*/
void delete_bin(std::vector<uint64_t> const & delete_files, estimate_ibf_arguments & ibf_args, std::filesystem::path path_in, bool samplewise);
