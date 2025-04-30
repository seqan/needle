// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "misc/get_expression_thresholds.hpp"

#include <cassert>

// Calculate expression thresholds and sizes
// Calculate expression thresholds by taking median recursively
void get_expression_thresholds(uint8_t const number_expression_thresholds,
                               robin_hood::unordered_node_map<uint64_t, uint16_t> const & hash_table,
                               std::vector<uint16_t> & expression_thresholds,
                               std::vector<uint64_t> & sizes,
                               robin_hood::unordered_set<uint64_t> const & genome,
                               uint8_t const cutoff,
                               bool const all)
{
    expression_thresholds.clear();
    sizes.clear();

    std::vector<uint16_t> counts;
    for (auto && elem : hash_table)
    {
        if (all || genome.contains(elem.first))
            counts.push_back(elem.second);
    }
    size_t const num_counts = counts.size();

    size_t median_fraction{2}; // First median at num_counts / 2, second one at num_counts / 2 + num_counts / 4, etc.
    size_t prev_median_position{};
    size_t median_position{};
    uint16_t prev_expression{};
    uint16_t expression{};
    uint16_t const max_count = std::ranges::max(counts);

    // Zero Level = cutoff + 1
    expression_thresholds.push_back(cutoff + 1u);

    while (expression_thresholds.size() < number_expression_thresholds && prev_expression < max_count
           && median_fraction < num_counts)
    {
        median_position = prev_median_position + num_counts / median_fraction;
        std::nth_element(counts.begin() + prev_median_position, counts.begin() + median_position, counts.end());
        expression = counts[median_position];
        prev_median_position = median_position;
        median_fraction = median_fraction * 2;

        // If expression does not change compared to previous one, do not store it again as an expression threshold.
        if (expression != prev_expression || expression_thresholds.empty())
        {
            expression_thresholds.push_back(expression);
            sizes.push_back(prev_median_position);
        }

        prev_expression = expression;
    }
    sizes.push_back(prev_expression);
    // In case not all levels have a threshold, give the last levels a maximal threshold, which can not be met by any minimiser.
    expression_thresholds.resize(number_expression_thresholds, max_count + 1u);

    assert(counts.size() == num_counts);
}
