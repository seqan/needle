// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "misc/read_levels.hpp"

#include <cassert>
#include <charconv>
#include <fstream>
#include <string>

#include "misc/charconv" // Can be dropped once clang-22 is released and clang-19 dropped.
#include "misc/needle_matrix.hpp"

// Reads the level file ibf creates
template <typename float_or_int>
void read_levels(std::vector<std::vector<float_or_int>> & expressions, std::filesystem::path const & filename)
{
    std::ifstream fin{filename};
    expressions.clear();

    std::string line;
    while (std::getline(fin, line))
    {
        // Strict format: The only line starting with '/' is the terminator.
        if (line.starts_with('/'))
            break;

        // After the last number, there is an extra space separator.
        if (line.ends_with(' '))
            line.pop_back();

        expressions.emplace_back();
        auto & row = expressions.back();

        char const * ptr = line.data();
        char const * end = ptr + line.size();

        while (ptr < end)
        {
            // Strict format: We are either at the start of a number or a single space separator.
            if (*ptr == ' ')
                ++ptr;

            float_or_int value{};
            auto res = std::from_chars(ptr, end, value);

            // Debug Mode: Verify parsing succeeded.
            // Release Mode: Assume success (res.ec == std::errc()).
            assert(res.ec == std::errc());

            row.push_back(value);
            ptr = res.ptr;
        }
    }
}

template void read_levels<uint16_t>(std::vector<std::vector<uint16_t>> & expressions,
                                    std::filesystem::path const & filename);
template void read_levels<double>(std::vector<std::vector<double>> & expressions,
                                  std::filesystem::path const & filename);

// Overload for needle_matrix
template <typename float_or_int>
void read_levels(needle_matrix<float_or_int> & expressions, std::filesystem::path const & filename)
{
    std::ifstream fin{filename};
    std::vector<float_or_int> tmp;
    size_t levels{};

    std::string line;
    while (std::getline(fin, line))
    {
        // Strict format: The only line starting with '/' is the terminator.
        if (line.starts_with('/'))
            break;

        // After the last number, there is an extra space separator.
        if (line.ends_with(' '))
            line.pop_back();

        char const * ptr = line.data();
        char const * end = ptr + line.size();

        while (ptr < end)
        {
            // Strict format: We are either at the start of a number or a single space separator.
            if (*ptr == ' ')
                ++ptr;

            float_or_int value{};
            auto res = std::from_chars(ptr, end, value);

            // Debug Mode: Verify parsing succeeded.
            // Release Mode: Assume success (res.ec == std::errc()).
            assert(res.ec == std::errc());

            tmp.push_back(value);
            ptr = res.ptr;
        }

        ++levels;
    }

    size_t experiments = tmp.size() / levels;
    assert(experiments * levels == tmp.size());
    expressions = needle_matrix<float_or_int>{std::move(tmp), levels, experiments};
}

template void read_levels<uint16_t>(needle_matrix<uint16_t> & expressions, std::filesystem::path const & filename);
template void read_levels<double>(needle_matrix<double> & expressions, std::filesystem::path const & filename);
