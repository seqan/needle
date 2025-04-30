// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "misc/read_levels.hpp"

#include <fstream>

#include <seqan3/io/views/detail/istreambuf_view.hpp>
#include <seqan3/io/views/detail/take_until_view.hpp>
#include <seqan3/utility/char_operations/predicate.hpp>

template <typename float_or_int>
float_or_int parse(std::string const & buffer)
{
    // std::from_chars would be better, but float version not available for every compiler
    if constexpr (std::same_as<uint16_t, float_or_int>)
        return std::stoi(buffer);
    else
        return std::stod(buffer);
}

// Reads the level file ibf creates
template <typename float_or_int>
void read_levels(std::vector<std::vector<float_or_int>> & expressions, std::filesystem::path filename)
{
    std::ifstream fin{filename};
    auto stream_view = seqan3::detail::istreambuf(fin);
    auto stream_it = std::ranges::begin(stream_view);
    size_t j{};
    std::vector<float_or_int> empty_vector{};
    std::string buffer{};

    // Read line = expression levels
    do
    {
        buffer.clear();

        if (j == expressions.size())
            expressions.push_back(empty_vector);

        std::ranges::copy(stream_view | seqan3::detail::take_until_or_throw(seqan3::is_char<' '>),
                          std::back_inserter(buffer));

        expressions[j].push_back(parse<float_or_int>(buffer));

        if (*stream_it != '/')
            ++stream_it;

        if (*stream_it == '\n')
        {
            ++stream_it;
            j++;
        }
    }
    while (*stream_it != '/');
}

template void read_levels<uint16_t>(std::vector<std::vector<uint16_t>> & expressions, std::filesystem::path filename);
template void read_levels<double>(std::vector<std::vector<double>> & expressions, std::filesystem::path filename);
