// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstdint>
#include <robin_hood.h>
#include <vector>

// Calculate expression thresholds and sizes
// Calculate expression thresholds by taking median recursively
void get_expression_thresholds(uint8_t const number_expression_thresholds,
                               robin_hood::unordered_node_map<uint64_t, uint16_t> const & hash_table,
                               std::vector<uint16_t> & expression_thresholds,
                               std::vector<uint64_t> & sizes,
                               robin_hood::unordered_set<uint64_t> const & genome,
                               uint8_t const cutoff,
                               bool const all = true);
