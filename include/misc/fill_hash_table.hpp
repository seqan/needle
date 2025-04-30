// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstdint>
#include <robin_hood.h>

#include "shared.hpp"

// Fill hash table with minimisers greater than the cutoff.
void fill_hash_table(minimiser_arguments const & args,
                     sequence_file_t & fin,
                     robin_hood::unordered_node_map<uint64_t, uint16_t> & hash_table,
                     robin_hood::unordered_node_map<uint64_t, uint8_t> & cutoff_table,
                     robin_hood::unordered_set<uint64_t> const & include_set_table,
                     robin_hood::unordered_set<uint64_t> const & exclude_set_table,
                     bool const only_include = false,
                     uint8_t cutoff = 0);

void fill_hash_table_parallel(minimiser_arguments const & args,
                              sequence_file_t & fin,
                              robin_hood::unordered_node_map<uint64_t, uint16_t> & hash_table,
                              robin_hood::unordered_node_map<uint64_t, uint8_t> & cutoff_table,
                              robin_hood::unordered_set<uint64_t> const & include_set_table,
                              robin_hood::unordered_set<uint64_t> const & exclude_set_table,
                              bool const only_include = false,
                              uint8_t cutoff = 0);
