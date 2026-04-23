// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstddef>
#include <ranges>
#include <span>
#include <vector>

template <typename value_t>
class needle_matrix
{
    std::vector<value_t> data_;
    size_t levels_{};
    size_t experiments_{};

public:
    needle_matrix() = default;
    needle_matrix(needle_matrix const &) = default;
    needle_matrix(needle_matrix &&) = default;
    needle_matrix & operator=(needle_matrix const &) = default;
    needle_matrix & operator=(needle_matrix &&) = default;
    ~needle_matrix() = default;

    needle_matrix(std::vector<value_t> data, size_t levels, size_t experiments) :
        data_(std::move(data)),
        levels_(levels),
        experiments_(experiments)
    {}

    needle_matrix(size_t levels, size_t experiments) :
        data_(levels * experiments),
        levels_(levels),
        experiments_(experiments)
    {}

    template <typename self_t>
    [[nodiscard]] constexpr auto * data(this self_t && self) noexcept
    {
        return self.data_.data();
    }

    // Access via [level, experiment]
    template <typename self_t>
    [[nodiscard]] constexpr auto && operator[](this self_t && self, size_t lvl, size_t exp) noexcept
    {
        return self[lvl * self.experiments() + exp];
    }

    // Flat 1D access (for bin-wise operations)
    template <typename self_t>
    [[nodiscard]] constexpr auto && operator[](this self_t && self, size_t bin) noexcept
    {
        return self.data()[bin];
    }

    // Returns a contiguous span for a single level (row)
    template <typename self_t>
    [[nodiscard]] constexpr auto level(this self_t && self, size_t lvl) noexcept
    {
        return std::span(self.data() + (lvl * self.experiments()), self.experiments());
    }

    // Returns a transform view for a single experiment (column)
    template <typename self_t>
    [[nodiscard]] constexpr auto experiment(this self_t && self, size_t exp) noexcept
    {
        return std::views::iota(size_t{0}, self.levels())
             | std::views::transform(
                   [&self, exp](size_t lvl) -> auto &
                   {
                       return self[lvl, exp];
                   });
    }

    constexpr void zero() noexcept
    {
        std::ranges::fill_n(data(), size(), value_t{});
    }

    [[nodiscard]] constexpr size_t levels() const noexcept
    {
        return levels_;
    }
    [[nodiscard]] constexpr size_t experiments() const noexcept
    {
        return experiments_;
    }
    [[nodiscard]] constexpr size_t size() const noexcept
    {
        return data_.size();
    }
};
