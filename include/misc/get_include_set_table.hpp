// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <vector>

#include <hibf/contrib/robin_hood.hpp>

#include "shared.hpp"

// Create set with hashes from the minimisers from an include or exclude file.
void get_include_set_table(minimiser_arguments const & args,
                           std::filesystem::path const include_file,
                           robin_hood::unordered_set<uint64_t> & include_table);
