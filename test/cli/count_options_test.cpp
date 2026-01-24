// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <string> // strings

#include "../app_test.hpp"

struct count_options_test : public app_test
{};

TEST_F(count_options_test, no_options)
{
    app_test_result result = execute_app("count");
    std::string expected{"needle-count - Get expression value depending on minimizers. This function is an alternative "
                         "to pseudoaligners like kallisto. It estimates the expression value for all sequences in the "
                         "genome file based on the exact minimiser occurrences of the given sequence files. Please run "
                         "genome beforehand to create the genome "
                         "file.\n======================================================================================"
                         "============================================================================================="
                         "============================================================================================="
                         "==================================================\n"
                         "    needle count [-p|--paired] [-k|--kmer uint8] [-w|--window uint32] [--shape\n"
                         "    uint64] [--seed uint64] [-o|--out path] [-t|--threads uint16] --include\n"
                         "    path [--genome path] [--] path...\n"
                         "    Try -h or --help for more information.\n"};
    EXPECT_SUCCESS(result);
#ifndef __APPLE__ // uint64_t vs unsigned long in sharg 1.2.1
    EXPECT_EQ(result.out, expected);
#endif
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(count_options_test, fail_no_argument)
{
    app_test_result result =
        execute_app("count", "--seed 0 --genome", data("mini_gen.genome"), " --include", data("mini_gen.fasta"));
    std::string expected{
        "[Error] Not enough positional arguments provided (Need at least 1). See -h/--help for more information.\n"};
    EXPECT_FAILURE(result);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(count_options_test, with_arguments)
{
    app_test_result result = execute_app("count -k 4 -w 4 --seed 0 --genome",
                                         data("mini_gen.genome"),
                                         " --include",
                                         data("mini_gen.fasta"),
                                         data("mini_example.fasta"));
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(count_options_test, multithreads)
{
    app_test_result result = execute_app("count -k 4 -w 8 -t 2 --genome",
                                         data("mini_gen.genome"),
                                         " --include",
                                         data("mini_gen.fasta"),
                                         data("mini_example.fasta"),
                                         data("mini_example.fasta"));
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(count_options_test, paired)
{
    app_test_result result = execute_app("count -k 4 -w 8 -p --genome",
                                         data("mini_gen.genome"),
                                         " --include",
                                         data("mini_gen.fasta"),
                                         data("mini_example.fasta"),
                                         data("mini_example.fasta"));
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}

struct genome_options_test : public app_test
{};

TEST_F(genome_options_test, no_options)
{
    app_test_result result = execute_app("genome");
    std::string expected{"needle-genome - Creates the genome file necessary as an input to count.\n"
                         "=======================================================================\n"
                         "    needle genome [-k|--kmer uint8] [-w|--window uint32] [--shape uint64]\n"
                         "    [--seed uint64] [-o|--out path] [-t|--threads uint16] [--exclude path]\n"
                         "    [--] path\n"
                         "    Try -h or --help for more information.\n"};
    EXPECT_SUCCESS(result);
#ifndef __APPLE__ // uint64_t vs unsigned long in sharg 1.2.1
    EXPECT_EQ(result.out, expected);
#endif
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(genome_options_test, fail_no_argument)
{
    app_test_result result = execute_app("genome", "--seed 0");
    std::string expected{
        "[Error] Not enough positional arguments provided (Need at least 1). See -h/--help for more information.\n"};
    EXPECT_FAILURE(result);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(genome_options_test, with_arguments)
{
    app_test_result result = execute_app("genome -k 4 -w 4 --seed 0 ", data("mini_gen.fasta"));
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}
