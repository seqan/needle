// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "misc/check_for_fasta_format.hpp"

#include <algorithm>

// Check if file has fasta format to estimate cutoffs.
bool check_for_fasta_format(std::vector<std::string> const & valid_extensions, std::string const & file_path)
{

    auto case_insensitive_string_ends_with = [&](std::string_view str, std::string_view suffix)
    {
        size_t const suffix_length{suffix.size()};
        size_t const str_length{str.size()};
        return suffix_length > str_length ? false
                                          : std::ranges::equal(str.substr(str_length - suffix_length),
                                                               suffix,
                                                               [](char const chr1, char const chr2)
                                                               {
                                                                   return std::tolower(chr1) == std::tolower(chr2);
                                                               });
    };

    auto case_insensitive_ends_with = [&](std::string const & ext)
    {
        return case_insensitive_string_ends_with(file_path, ext);
    };

    return std::ranges::find_if(valid_extensions, case_insensitive_ends_with) != valid_extensions.end();
}
