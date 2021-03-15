#pragma once

#include <seqan3/std/filesystem>

#include "minimiser.h"

/*!\brief The arguments necessary for a search.
 * \param std::filesystem::path search_file The sequence file containing the transcripts to be searched for.
 * \param std::filesystem::path exp_file    The expression file (a tab seperated file containing expression information
 *                                          per transcript given). Optimally, it is ordered by the expression levels.
 *                                          Does not need to be specified.
 * \param std::filesystem::path path_in     The path to the directory where the IBFs can be found. Default: Current
 *                                          directory.
 * \param float expression                  The expression level that should be used when searching for a transcript
 *                                          (if no expression file is given to specify this individually for
 *                                          different transcripts).
 * \param float threshold                   The percentage of how many minimisers of a transcript need to be found in an
 *                                          IBF to be considered as contained in a certain IBF.
 *
 */
struct estimate_arguments
{
    std::filesystem::path search_file;
    std::filesystem::path path_in{"./"};
    std::vector<uint32_t> expressions{};
    float threshold{0.5};

};

/*! \brief Function to estimate expression value.
*  \param args        The arguments.
*  \param estimate_args The search arguments.
*  \param ibf         The ibf determing what kind ibf is used (compressed or uncompressed).
*  \param file_out    The file where results should be stored to.
*  \param search_file The sequence file with the sequences which expression value should be estimated.
*  \param path_in     The directory where the ibfs can be found.
*/
template <class IBFType>
void estimate(arguments const & args, estimate_arguments const & estimate_args, IBFType & ibf, std::filesystem::path file_out,
              std::filesystem::path search_file, std::filesystem::path path_in);

/*! \brief Function to estimate expression value.
*  \param args        The arguments.
*  \param estimate_args The search arguments.
*  \param ibf         The ibf determing what kind ibf is used (compressed or uncompressed).
*  \param file_out    The file where results should be stored to.
*  \param search_file The sequence file with the sequences which expression value should be estimated.
*  \param path_in     The directory where the ibfs can be found.
*  \param level_file Path to the header files, where expression levels can be found.
*/
template <class IBFType>
void estimate(arguments const & args, estimate_arguments const & estimate_args, IBFType & ibf, std::filesystem::path file_out,
            std::filesystem::path search_file, std::filesystem::path path_in, std::filesystem::path level_file);

/*! \brief Function, which calls the estimate function.
*  \param args        The arguments.
*  \param estimate_args The search arguments.
*  \param file_out    The file where results should be stored to.
*  \param search_file The sequence file with the sequences which expression value should be estimated.
*  \param path_in     The directory where the ibfs can be found.
*  \param level_file  Path to the header files, where expression levels can be found.
*/
void call_estimate(arguments const & args, estimate_arguments const & estimate_args, std::filesystem::path file_out,
            std::filesystem::path search_file, std::filesystem::path path_in,std::filesystem::path level_file);
