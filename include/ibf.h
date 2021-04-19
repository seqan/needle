#pragma once

#include <math.h>
#include <numeric>
#include <string>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/container/concatenated_sequences.hpp>
#include <seqan3/std/filesystem>

#include "minimiser.h"

//!\brief specific arguments needed for constructing an IBF
struct ibf_arguments
{
    std::vector<std::filesystem::path> sequence_files;
    std::filesystem::path include_file; // Needs to be defined when only minimizers appearing in this file should be stored
    std::filesystem::path exclude_file; // Needs to be defined when minimizers appearing in this file should NOT be stored
    std::vector<size_t> bin_size{}; // The bin size of one IBF, can be different for different expression levels
    size_t num_hash{1}; // Number of hash functions to use, default 1
    std::filesystem::path path_out{"./"}; // Path where IBFs should be stored
    std::vector<uint32_t> expression_levels{}; // Expression levels which should be created
    std::vector<int> samples{}; // Can be used to indicate that sequence files belong to the same experiment
    bool paired = false; // If true, than experiments are seen as paired-end experiments
    // Which expression values should be ignored during calculation of the normalization_method, default is zero
    std::vector<uint32_t> cutoffs{};
    //std::string expression_method{"median"}; // Method to calculate expression levels
    bool experiment_names = false; // Flag, if names of experiment should be stored in a txt file
    uint8_t number_expression_levels{};
    // Flag, if true the "median" approach is used to determine the expression_levels individually for a sample
    bool set_expression_levels_samplewise{false};
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
void get_minimisers(arguments const & args, seqan3::concatenated_sequences<seqan3::dna4_vector> const & sequences,
                    robin_hood::unordered_node_map<uint64_t, uint16_t> & hash_table,
                    robin_hood::unordered_set<uint64_t> const & genome_set_table,
                    std::filesystem::path const & genome_file = "", bool only_genome = false);

/*!\brief Get the concrete expression values (= median of all counts of one transcript).
* \param args               The minimiser arguments to use (seed, shape, window size).
* \param sequence_files     The sequence files, which contains the reads.
* \param genome_file        A file containing the transcripts which expression values should be determined.
* \param out_path        The output path, where results are stored
* \param paired             Flag to indicate if input data is paired or not.
*/
void count(arguments const & args, std::vector<std::filesystem::path> sequence_files, std::filesystem::path genome_file,
           std::filesystem::path out_path, bool paired);

/*!\brief Set arguments for creating IBF.
* \param args               The minimiser arguments to use (seed, shape, window size).
* \param ibf_args           The IBF specific arguments to use (bin size, number of hash functions, ...). See
*                           struct ibf_arguments.
* \param genome_sequences   Data structure, where the sequences of the genome file are stored.
* \param genome_set_table   Data structure, where the minimisers found in a genome mask are stored.
*/
void set_arguments_ibf(arguments const & args, ibf_arguments & ibf_args,
                   seqan3::concatenated_sequences<seqan3::dna4_vector> & genome_sequences,
                   robin_hood::unordered_set<uint64_t> & genome_set_table);

/*!\brief Set arguments for creating IBF.
* \param ibf_args           The IBF specific arguments to use (bin size, number of hash functions, ...). See
*                           struct ibf_arguments.
*/
void set_arguments(ibf_arguments & ibf_args);

/*!\brief Reads a binary file function minimiser creates
* \param hash_table         The hash table to store minimisers into.
* \param filename           The filename of the binary file.
*/
void read_binary(robin_hood::unordered_node_map<uint64_t, uint16_t> & hash_table, std::filesystem::path filename);

/*!\brief Reads a header file function minimiser creates
* \param args                 The minimiser arguments to use (seed, shape, window size).
* \param ibf_args             The IBF specific arguments to use (bin size, number of hash functions, ...). See
*                             struct ibf_arguments.
* \param filename             The filename of the binary file.
* \param counts               Vector, where the number of minimiser at a certain expression level should be stored into.
*/
void read_header(arguments & args, ibf_arguments & ibf_args, std::filesystem::path filename,
                 std::vector<uint16_t> & counts);

/*! \brief Create IBF.
 * \param args         The minimiser arguments to use (seed, shape, window size).
 * \param ibf_args     The IBF specific arguments to use (bin size, number of hash functions, ...). See
 *                     struct ibf_arguments.
 *  \returns The normalized expression values per experiment.
 */
std::vector<uint32_t> ibf(arguments const & args, ibf_arguments & ibf_args);

/*! \brief Create IBF based on the minimiser and header files
 * \param minimiser_files  A vector of minimiser file paths.
 * \param args             The minimiser arguments to use (seed, shape, window size).
 * \param ibf_args         The IBF specific arguments to use (bin size, number of hash functions, ...). See
 *                         struct ibf_arguments.
 *  \returns The normalized expression values per experiment.
 */
std::vector<uint32_t> ibf(std::vector<std::filesystem::path> minimiser_files, arguments & args,
                          ibf_arguments & ibf_args);

void minimiser(arguments const & args, ibf_arguments & ibf_args);

void build_ibf(arguments & args, ibf_arguments & ibf_args, float fpr = 0.05);
