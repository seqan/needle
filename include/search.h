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
struct search_arguments
{
    std::filesystem::path search_file;
    std::filesystem::path exp_file;
    std::filesystem::path path_in{"./"};
    float expression{1.0};
    float threshold{0.5};

};

/*! \brief Function, which searches for transcripts in IBF of a given expression level.
*  \param ibf         The IBF.
*  \param args        The arguments.
*  \param search_args The search arguments.
* \returns result vector.
*/
template <class IBFType>
std::vector<uint32_t> do_search(IBFType & ibf, arguments const & args, search_arguments const & search_args);

/*! \brief Function, which calls the search functions.
*  \param args        The arguments.
*  \param search_args The search arguments.
* \returns result vector.
*/
std::vector<uint32_t> search(arguments const & args, search_arguments const & search_args);
