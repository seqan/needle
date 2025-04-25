// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <string> // strings

#include "../app_test.hpp"

struct needle_options_test : public app_test
{};

TEST_F(needle_options_test, no_options)
{
    app_test_result result = execute_app();
    std::string expected{"needle\n"
                         "======\n"
                         "    Try -h or --help for more information.\n"};
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(needle_options_test, fail_no_argument)
{
    app_test_result result = execute_app("-v");
    std::string expected{"[Error] Unknown option -v. In case this is meant to be a non-option/argument/parameter, "
                         "please specify the start of non-options with '--'. See -h/--help for program information.\n"};
    EXPECT_FAILURE(result);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}
