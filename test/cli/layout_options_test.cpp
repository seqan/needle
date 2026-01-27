// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <string> // strings

#include "../app_test.hpp"

struct needle_layout_test : public app_test
{};

TEST_F(needle_layout_test, no_options)
{
    app_test_result result = execute_app("layout");
    std::string expected{"needle-layout - Compute an HIBF layout\n"
                         "======================================\n"
                         "     --input <file> [--output <file>] [--threads <number>] [--kmer <number>]\n"
                         "    [--fpr <number>] [--hash <number>] [--disable-estimate-union]\n"
                         "    [--disable-rearrangement]\n"
                         "    Try -h or --help for more information.\n"};
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}
