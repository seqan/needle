// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <string>
#include <vector>

// Check if file has fasta format to estimate cutoffs.
bool check_for_fasta_format(std::vector<std::string> const & valid_extensions, std::string const & file_path);
