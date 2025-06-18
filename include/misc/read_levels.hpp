// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <filesystem>
#include <vector>

// Reads the level file ibf creates
template <typename float_or_int>
void read_levels(std::vector<std::vector<float_or_int>> & expressions, std::filesystem::path const & filename);
