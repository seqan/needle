// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <vector>

// Check and set samples and cutoffs
void check_cutoffs_samples(std::vector<std::filesystem::path> const & sequence_files,
                           bool const paired,
                           std::vector<size_t> & samples,
                           std::vector<uint8_t> & cutoffs);
