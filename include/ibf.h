#pragma once

#include <iostream>
#include <math.h>
#include <numeric>
#include <string>

#include <seqan3/alphabet/container/concatenated_sequences.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/std/filesystem>

#include "shared.h"

struct minimiser_arguments
{
    std::filesystem::path include_file; // Needs to be defined when only minimizers appearing in this file should be stored
    std::filesystem::path exclude_file; // Needs to be defined when minimizers appearing in this file should NOT be stored
    std::vector<int> samples{}; // Can be used to indicate that sequence files belong to the same experiment
    bool paired = false; // If true, than experiments are seen as paired-end experiments
    std::vector<uint8_t> cutoffs{};
    bool experiment_names = false; // Flag, if names of experiment should be stored in a txt file
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

/*!\brief Calculate best bin size based on number of elements maximal inserted, false positive rate and number of
 *        hash functions. See: https://hur.st/bloomfilter/
 * \param count     The number of elements to be stored.
 * \param fpr       The false positive rate to use.
 * \param num_hash  The number of hash functions to use.
 * \returns bin_size
 */
inline uint64_t get_bin_size(uint64_t count, float fpr, size_t num_hash)
{
    return std::ceil((count * std::log(fpr)) / std::log(1 / std::pow(2, std::log(2))));
}

/*!\brief Gets all sequences from a specified number of sequence files.
 * \param sequence_files A vector of paths to the sequence files.
 * \param sequences      The data strucuture, where the sequenecs should be stored in.
 * \param min_len        The minimum length a sequence should have to be stored.
 * \param first          The first position of in the sequence file vector, which should be considered. Default: 0.
 * \param num_exp        The number of sequence files that should be considered. Default: 1.
 */
void get_sequences(std::vector<std::filesystem::path> const & sequence_files,
                   seqan3::concatenated_sequences<seqan3::dna4_vector> & sequences, uint16_t min_len,
                   unsigned first = 0, unsigned num_exp = 1);

/*!\brief Gets all minimisers from all sequences in a given vector.
* \param args                      The minimiser arguments to use (seed, shape, window size).
* \param sequences                 The data strucuture, where the sequenecs are stored.
* \param hash_table                The hash table, where minimisers should be stored.
* \param genome_set_table          The minimisers found in a genome mask.
* \param genome_file               The file to the genome mask. Default: "".
* \param only_genome               True, if only minimisers found in the genome mask should be stored. Default: False.
*/
void get_minimisers(min_arguments const & args, seqan3::concatenated_sequences<seqan3::dna4_vector> const & sequences,
                    robin_hood::unordered_node_map<uint64_t, uint16_t> & hash_table,
                    robin_hood::unordered_set<uint64_t> const & genome_set_table,
                    std::filesystem::path const & genome_file = "", bool only_genome = false);

/*!\brief Get the concrete expression values (= median of all counts of one transcript) for given experiments.
*         This function can be used to estimate how good the median approach can be, if all count values are available.
* \param args               The minimiser arguments to use (seed, shape, window size).
* \param sequence_files     The sequence files, which contains the reads.
* \param genome_file        A file containing the transcripts which expression values should be determined.
* \param exclude_file       A file containing minimizers which should be ignored.
* \param paired             Flag to indicate if input data is paired or not.
*/
void count(min_arguments const & args, std::vector<std::filesystem::path> sequence_files, std::filesystem::path genome_file,
           std::filesystem::path exclude_file, bool paired);

/*!\brief Reads a binary file that needle minimiser creates.
* \param filename           The filename of the binary file.
* \param hash_table         The hash table to store minimisers into.

*/
void read_binary(std::filesystem::path filename, robin_hood::unordered_node_map<uint64_t, uint16_t> & hash_table);

/*!\brief Reads the beginning of a binary file that needle minimiser creates.
* \param args               Min arguments.
* \param filename           The filename of the binary file.
* \param num_of_minimisers  Variable, where to number of minimisers should be stored.

*/
void read_binary_start(min_arguments & args, std::filesystem::path filename, uint64_t & num_of_minimisers);

/*! \brief Create IBF.
 * \param sequence_files  A vector of sequence file paths.
 * \param ibf_args        The IBF specific arguments to use (bin size, number of hash functions, ...). See
 *                        struct ibf_arguments.
 * \param minimiser_args  The minimiser specific arguments to use.
 * \param fpr             The average false positive rate that should be used.
 * \param expression_by_genome_file File that contains the only minimisers that should be comnsidered for the
 *                                  determination of the expression_levels.
 * \param num_hash        The number of hash functions to use.
 *  \returns The normalized expression values per experiment.
 */
std::vector<uint16_t> ibf(std::vector<std::filesystem::path> const & sequence_files, estimate_ibf_arguments & ibf_args,
                          minimiser_arguments & minimiser_args, std::vector<double> & fpr,
                          std::filesystem::path const expression_by_genome_file = "",
                          size_t num_hash = 1);

/*! \brief Create IBF based on the minimiser and header files
 * \param minimiser_files  A vector of minimiser file paths.
 * \param ibf_args         The IBF specific arguments to use (bin size, number of hash functions, ...). See
 *                         struct ibf_arguments.
 * \param fpr             The average false positive rate that should be used.
 * \param expression_by_genome_file File that contains the only minimisers that should be comnsidered for the
 *                                  determination of the expression_levels.
 * \param num_hash        The number of hash functions to use.
 *  \returns The normalized expression values per experiment.
 */
std::vector<uint16_t> ibf(std::vector<std::filesystem::path> const & minimiser_files,
                          estimate_ibf_arguments & ibf_args, std::vector<double> & fpr,
                          std::filesystem::path const expression_by_genome_file = "",
                          size_t num_hash = 1);

/*! \brief Create minimiser and header files.
* \param sequence_files  A vector of sequence file paths.
* \param args             The minimiser arguments to use (seed, shape, window size).
* \param minimiser_args  The minimiser specific arguments to use.
*/
void minimiser(std::vector<std::filesystem::path> const & sequence_files, min_arguments const & args,
               minimiser_arguments & minimiser_args);
