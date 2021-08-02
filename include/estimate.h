#pragma once

#include <seqan3/std/filesystem>

#include "shared.h"

/*!\brief The arguments necessary for a search.
 * \param std::filesystem::path search_file The sequence file containing the transcripts to be searched for.
 * \param std::filesystem::path path_in     The path to the directory where the IBFs can be found. Default: Current
 *                                          directory.
 * \param bool normalization_method         Flag, true if normalization should be used.
 *
 */
struct estimate_arguments
{
    std::filesystem::path search_file;
    std::filesystem::path path_in{"./"};
    // false: no normalization method, true: division by first expression value
    bool normalization_method{0};
};

/*! \brief Function, which calls the estimate function.
*  \param args          The arguments estimate and ibf use.
*  \param estimate_args The estimate arguments.
*  \param level_file    Path to the header files, where expression levels can be found.
*/
void call_estimate(estimate_ibf_arguments & args, estimate_arguments & estimate_args,
                   std::filesystem::path level_file = "");
