// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <iostream>

#include <seqan3/test/expect_range_eq.hpp>

#include "../app_test.hpp"
#include "estimate.hpp"
#include "ibf.hpp"
#include "shared.hpp"

// To prevent issues when running multiple CLI tests in parallel, give each CLI test unique names:
struct estimate_test : public app_test
{
    void initialization_args(estimate_ibf_arguments & args)
    {
        args.compressed = true;
        args.k = 4;
        args.shape = seqan3::ungapped{args.k};
        args.w_size = seqan3::window_size{4};
        args.s = seqan3::seed{0};
        args.path_out = "Estimate_Test_";
    }
};

TEST_F(estimate_test, small_example)
{
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    estimate_arguments estimate_args{};
    initialization_args(ibf_args);
    ibf_args.expression_thresholds = {1, 2, 4};
    std::vector<double> fpr = {0.05};
    std::vector<std::filesystem::path> sequence_files = {data("mini_example.fasta")};
    estimate_args.search_file = data("mini_gen.fasta");
    estimate_args.path_in = ibf_args.path_out;
    std::vector<uint8_t> cutoffs{};

    ibf(sequence_files, ibf_args, minimiser_args, fpr, cutoffs);
    ibf_args.path_out = "expression.out";
    call_estimate(ibf_args, estimate_args);

    std::ifstream output_file("expression.out");
    std::string line;
    std::string expected{"gen1\t3\t"};
    if (output_file.is_open())
    {
        while (std::getline(output_file, line))
        {
            EXPECT_EQ(expected, line);
        }
        output_file.close();
    }
}

TEST_F(estimate_test, small_example_uncompressed)
{
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    estimate_arguments estimate_args{};
    initialization_args(ibf_args);
    ibf_args.compressed = false;
    ibf_args.expression_thresholds = {1, 2, 4};
    std::vector<double> fpr = {0.05};
    std::vector<std::filesystem::path> sequence_files = {data("mini_example.fasta")};
    estimate_args.search_file = data("mini_gen.fasta");
    estimate_args.path_in = ibf_args.path_out;
    std::vector<uint8_t> cutoffs{};

    ibf(sequence_files, ibf_args, minimiser_args, fpr, cutoffs);
    ibf_args.path_out = "expression.out";
    call_estimate(ibf_args, estimate_args);

    std::ifstream output_file("expression.out");
    std::string line;
    std::string expected{"gen1\t3\t"};
    if (output_file.is_open())
    {
        while (std::getline(output_file, line))
        {
            EXPECT_EQ(expected, line);
        }
        output_file.close();
    }
}

TEST_F(estimate_test, small_example_gene_not_found)
{
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    estimate_arguments estimate_args{};
    initialization_args(ibf_args);
    ibf_args.expression_thresholds = {2, 4};
    estimate_args.search_file = data("mini_gen2.fasta");
    estimate_args.path_in = ibf_args.path_out;
    std::vector<double> fpr = {0.05};
    std::vector<std::filesystem::path> sequence_files = {data("mini_example.fasta")};
    std::vector<uint8_t> cutoffs{};

    ibf(sequence_files, ibf_args, minimiser_args, fpr, cutoffs);
    ibf_args.path_out = "expression.out";
    call_estimate(ibf_args, estimate_args);

    std::ifstream output_file("expression.out");
    std::string line;
    std::string expected{"gen2\t0\t"};
    if (output_file.is_open())
    {
        while (std::getline(output_file, line))
        {
            EXPECT_EQ(expected, line);
        }
        output_file.close();
    }
}

TEST_F(estimate_test, small_example_different_expressions_per_level)
{
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    estimate_arguments estimate_args{};
    initialization_args(ibf_args);
    ibf_args.number_expression_thresholds = 2;
    std::vector<double> fpr = {0.05};
    std::vector<uint8_t> cutoffs{};
    std::vector<std::filesystem::path> sequence_files = {data("mini_example.fasta")};

    minimiser(sequence_files, ibf_args, minimiser_args, cutoffs);
    std::vector<std::filesystem::path> minimiser_files{"Estimate_Test_mini_example.minimiser"};
    ASSERT_TRUE(std::filesystem::exists(minimiser_files[0]));
    ibf_args.expression_thresholds = {};
    ibf(minimiser_files, ibf_args, fpr);

    ibf_args.expression_thresholds = {0, 1, 2};
    estimate_args.search_file = data("mini_gen.fasta");
    estimate_args.path_in = ibf_args.path_out;
    ibf_args.path_out = "Test0_expression.out";
    call_estimate(ibf_args, estimate_args);

    std::ifstream output_file("Test0_expression.out");
    std::string line;
    std::string expected{"gen1\t3\t"};
    if (output_file.is_open())
    {
        while (std::getline(output_file, line))
        {
            EXPECT_EQ(expected, line);
        }
        output_file.close();
    }
}

TEST_F(estimate_test, small_example_different_expressions_per_level_normalization_1)
{
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    estimate_arguments estimate_args{};
    estimate_args.normalization_method = 1;
    initialization_args(ibf_args);
    ibf_args.number_expression_thresholds = 2;
    std::vector<double> fpr = {0.05};
    std::vector<uint8_t> cutoffs{};
    std::vector<std::filesystem::path> sequence_files = {data("mini_example.fasta")};

    minimiser(sequence_files, ibf_args, minimiser_args, cutoffs);
    std::vector<std::filesystem::path> minimiser_files{"Estimate_Test_mini_example.minimiser"};
    ASSERT_TRUE(std::filesystem::exists(minimiser_files[0]));
    ibf_args.expression_thresholds = {};
    ibf(minimiser_files, ibf_args, fpr);

    ibf_args.expression_thresholds = {0, 1, 2};
    estimate_args.search_file = data("mini_gen.fasta");
    estimate_args.path_in = ibf_args.path_out;
    ibf_args.path_out = "expression.out";
    call_estimate(ibf_args, estimate_args);

    std::ifstream output_file("expression.out");
    std::string line;
    std::string expected{"gen1\t1\t"};
    if (output_file.is_open())
    {
        while (std::getline(output_file, line))
        {
            EXPECT_EQ(expected, line);
        }
        output_file.close();
    }
}

TEST_F(estimate_test, small_example_different_expressions_per_level_normalization_1_uncompressed)
{
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    estimate_arguments estimate_args{};
    estimate_args.normalization_method = 1;
    initialization_args(ibf_args);
    ibf_args.number_expression_thresholds = 2;
    ibf_args.compressed = false;
    std::vector<double> fpr = {0.05};
    std::vector<uint8_t> cutoffs{};
    std::vector<std::filesystem::path> sequence_files = {data("mini_example.fasta")};

    minimiser(sequence_files, ibf_args, minimiser_args, cutoffs);
    std::vector<std::filesystem::path> minimiser_files{"Estimate_Test_mini_example.minimiser"};
    ASSERT_TRUE(std::filesystem::exists(minimiser_files[0]));
    ibf_args.expression_thresholds = {};
    ibf(minimiser_files, ibf_args, fpr);

    ibf_args.expression_thresholds = {0, 1, 2};
    estimate_args.search_file = data("mini_gen.fasta");
    estimate_args.path_in = ibf_args.path_out;
    ibf_args.path_out = "Test2_expression.out";
    call_estimate(ibf_args, estimate_args);

    std::ifstream output_file("Test2_expression.out");
    std::string line;
    std::string expected{"gen1\t1\t"};
    if (output_file.is_open())
    {
        while (std::getline(output_file, line))
        {
            EXPECT_EQ(expected, line);
        }
        output_file.close();
    }
}

TEST_F(estimate_test, example)
{
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    estimate_arguments estimate_args{};
    std::vector<double> fpr = {0.05};
    std::vector<std::filesystem::path> sequence_files = {data("exp_01.fasta"),
                                                         data("exp_02.fasta"),
                                                         data("exp_11.fasta"),
                                                         data("exp_12.fasta")};
    minimiser_args.samples = {2, 2};
    ibf_args.expression_thresholds = {4, 32};
    ibf_args.compressed = false;
    std::vector<uint8_t> cutoffs{};
    ibf(sequence_files, ibf_args, minimiser_args, fpr, cutoffs);

    estimate_args.search_file = data("gene.fasta");
    estimate_args.path_in = ibf_args.path_out;
    ibf_args.path_out = "Single_expression.out";
    call_estimate(ibf_args, estimate_args);

    std::ifstream output_file("Single_expression.out");
    std::string line;
    std::string expected{"GeneA\t10\t32\t"};
    if (output_file.is_open())
    {
        while (std::getline(output_file, line))
        {
            EXPECT_EQ(expected, line);
        }
        output_file.close();
    }
}

TEST_F(estimate_test, example_multiple_threads)
{
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    estimate_arguments estimate_args{};
    std::vector<std::filesystem::path> sequence_files = {data("exp_01.fasta"),
                                                         data("exp_02.fasta"),
                                                         data("exp_11.fasta"),
                                                         data("exp_12.fasta")};
    minimiser_args.samples = {2, 2};
    ibf_args.expression_thresholds = {4, 32};
    std::vector<double> fpr = {0.05};
    ibf_args.compressed = false;
    std::vector<uint8_t> cutoffs{};
    ibf(sequence_files, ibf_args, minimiser_args, fpr, cutoffs);
    ibf_args.threads = 2;

    estimate_args.search_file = data("gene4.fasta");
    estimate_args.path_in = ibf_args.path_out;
    ibf_args.path_out = "Multiple_expression.out";
    call_estimate(ibf_args, estimate_args);

    std::ifstream output_file("Multiple_expression.out");
    std::string line;
    std::string expected{"GeneA\t10\t32\t"};
    if (output_file.is_open())
    {
        while (std::getline(output_file, line))
        {
            EXPECT_EQ(expected, line);
        }
        output_file.close();
    }
}

TEST_F(estimate_test, example_different_expressions_per_level)
{
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    estimate_arguments estimate_args{};
    std::vector<std::filesystem::path> sequence_files = {data("exp_01.fasta"),
                                                         data("exp_02.fasta"),
                                                         data("exp_11.fasta"),
                                                         data("exp_12.fasta")};
    std::vector<uint8_t> cutoffs = {0, 0};
    minimiser_args.samples = {2, 2};
    ibf_args.number_expression_thresholds = 4;
    std::vector<double> fpr = {0.05};
    ibf_args.compressed = false;
    minimiser(sequence_files, ibf_args, minimiser_args, cutoffs);
    std::vector<std::filesystem::path> minimiser_files{"exp_01.minimiser", "exp_11.minimiser"};
    ASSERT_TRUE(std::filesystem::exists(minimiser_files[0]));
    ASSERT_TRUE(std::filesystem::exists(minimiser_files[1]));
    ibf(minimiser_files, ibf_args, fpr);

    ibf_args.expression_thresholds = {0, 1, 2};
    estimate_args.search_file = data("gene.fasta");
    estimate_args.path_in = ibf_args.path_out;
    ibf_args.path_out = "expression.out";
    call_estimate(ibf_args, estimate_args);

    std::ifstream output_file("expression.out");
    std::string line;
    // Count would expect 6 and 34
    std::string expected{"GeneA\t7\t26\t"};
    if (output_file.is_open())
    {
        while (std::getline(output_file, line))
        {
            EXPECT_EQ(expected, line);
        }
        output_file.close();
    }
}

TEST_F(estimate_test, example_different_expressions_per_level_multiple_threads)
{
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    estimate_arguments estimate_args{};
    std::vector<std::filesystem::path> sequence_files = {data("exp_01.fasta"),
                                                         data("exp_02.fasta"),
                                                         data("exp_11.fasta"),
                                                         data("exp_12.fasta")};
    std::vector<uint8_t> cutoffs = {0, 0};
    minimiser_args.samples = {2, 2};
    ibf_args.number_expression_thresholds = 4;
    std::vector<double> fpr = {0.05};
    ibf_args.compressed = false;
    minimiser(sequence_files, ibf_args, minimiser_args, cutoffs);
    std::vector<std::filesystem::path> minimiser_files{"exp_01.minimiser", "exp_11.minimiser"};
    ASSERT_TRUE(std::filesystem::exists(minimiser_files[0]));
    ASSERT_TRUE(std::filesystem::exists(minimiser_files[1]));
    ibf_args.expression_thresholds = {};
    ibf(minimiser_files, ibf_args, fpr);

    ibf_args.threads = 2;
    ibf_args.expression_thresholds = {0, 1, 2};
    estimate_args.search_file = data("gene4.fasta");
    estimate_args.path_in = ibf_args.path_out;
    ibf_args.path_out = "expression.out";
    call_estimate(ibf_args, estimate_args);

    std::ifstream output_file("expression.out");
    std::string line;
    // Count would expect 6 and 34
    std::string expected{"GeneA\t7\t26\t"};
    if (output_file.is_open())
    {
        while (std::getline(output_file, line))
        {
            EXPECT_EQ(expected, line);
        }
        output_file.close();
    }
}
