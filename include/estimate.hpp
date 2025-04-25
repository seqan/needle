// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <filesystem>

#include "shared.hpp"

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
*/
void call_estimate(estimate_ibf_arguments & args, estimate_arguments & estimate_args);
