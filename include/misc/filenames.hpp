// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <filesystem>

#include "shared.hpp"

namespace filenames
{

constexpr std::string_view ibf_prefix = "IBF_";
constexpr std::string_view ibf_samplewise_prefix = "IBF_Level_";
static_assert(ibf_samplewise_prefix.starts_with(ibf_prefix));

constexpr std::string_view ibf_levels_name = "IBF_Levels.levels";
static_assert(ibf_levels_name.starts_with(ibf_prefix));

constexpr std::string_view ibf_fprs_name = "IBF_FPRs.fprs";
static_assert(ibf_fprs_name.starts_with(ibf_prefix));

constexpr std::string_view ibf_deleted_name = "IBF_Deleted";
static_assert(ibf_deleted_name.starts_with(ibf_prefix));

constexpr std::string_view ibf_data_name = "IBF_Data";
static_assert(ibf_data_name.starts_with(ibf_prefix));

constexpr std::string_view stored_files_name = "Stored_Files.txt";

// Note: base_path is taken by value.
//       We need a local copy for concat/+= anyway, see https://en.cppreference.com/w/cpp/filesystem/path/concat.
// We also want to avoid unnecessary copies of the path: https://godbolt.org/z/Mxrjxf8h3

[[nodiscard]] inline std::filesystem::path
ibf(std::filesystem::path base_path, bool const samplewise, uint16_t const level, estimate_ibf_arguments const & args)
{
    if (samplewise)
        base_path += (ibf_samplewise_prefix.data() + std::to_string(level));
    else
        base_path += (ibf_prefix.data() + std::to_string(args.expression_thresholds[level]));
    return base_path;
}

[[nodiscard]] inline std::filesystem::path levels(std::filesystem::path base_path)
{
    base_path += ibf_levels_name;
    return base_path;
}

[[nodiscard]] inline std::filesystem::path fprs(std::filesystem::path base_path)
{
    base_path += ibf_fprs_name;
    return base_path;
}

[[nodiscard]] inline std::filesystem::path data(std::filesystem::path base_path)
{
    base_path += ibf_data_name;
    return base_path;
}

[[nodiscard]] inline std::filesystem::path deleted(std::filesystem::path base_path)
{
    base_path += ibf_deleted_name;
    return base_path;
}

[[nodiscard]] inline std::filesystem::path stored(std::filesystem::path base_path)
{
    base_path += stored_files_name;
    return base_path;
}

} // namespace filenames
