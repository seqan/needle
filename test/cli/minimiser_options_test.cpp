// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <string> // strings

#include "../app_test.hpp"

struct minimiser_options_test : public app_test
{};

TEST_F(minimiser_options_test, no_options)
{
    app_test_result result = execute_app("minimiser");
    std::string expected{"needle-minimiser - Calculates minimiser for given experiments.\n"
                         "==============================================================\n"
                         "    Try -h or --help for more information.\n"};
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(minimiser_options_test, fail_no_argument)
{
    app_test_result result = execute_app("minimiser", "--seed 0");
    std::string expected{
        "[Error] Not enough positional arguments provided (Need at least 1). See -h/--help for more information.\n"};
    EXPECT_FAILURE(result);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(minimiser_options_test, with_arguments)
{
    app_test_result result = execute_app("minimiser -k 4 -w 4 --seed 0", data("mini_example.fasta"));
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(minimiser_options_test, cutoff)
{
    app_test_result result = execute_app("minimiser -k 4 -w 8 --cutoff 2", data("mini_example.fasta"));
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(minimiser_options_test, multiple_sample)
{
    app_test_result result =
        execute_app("minimiser -k 4 -w 8 --samples 2 ", data("mini_example.fasta"), data("mini_example.fasta"));
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(minimiser_options_test, multithreads)
{
    app_test_result result =
        execute_app("minimiser -k 4 -w 8 -t 2", data("mini_example.fasta"), data("mini_example.fasta"));
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(minimiser_options_test, paired)
{
    app_test_result result =
        execute_app("minimiser -k 4 -w 8 -p ", data("mini_example.fasta"), data("mini_example.fasta"));
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(minimiser_options_test, invalid_argument)
{
    app_test_result result =
        execute_app("minimiser -k 4 -w 8 --samples 3 ", data("mini_example.fasta"), data("mini_example.fasta"));
    std::string expected{"[Error] Incorrect command line input for multiple-samples.\n"};
    EXPECT_FAILURE(result);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, expected);
}
