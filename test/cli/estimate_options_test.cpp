// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <string> // strings

#include "../app_test.hpp"
#include "ibf.hpp"
#include "shared.hpp"

struct estimate_options_test : public app_test
{};

TEST_F(estimate_options_test, no_options)
{
    app_test_result result = execute_app("estimate");
    std::string expected{"needle-estimate - Estimate expression value of transcript based on the Needle index.\n"
                         "====================================================================================\n"
                         "    Try -h or --help for more information.\n"};
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(estimate_options_test, fail_no_argument)
{
    app_test_result result = execute_app("estimate", "-m");
    std::string expected{
        "[Error] Not enough positional arguments provided (Need at least 1). See -h/--help for more information.\n"};
    EXPECT_FAILURE(result);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(estimate_options_test, with_argument)
{
    estimate_ibf_arguments ibf_args{};
    minimiser_file_input_arguments minimiser_args{};
    ibf_args.expression_thresholds = {1, 2};
    std::vector<double> fpr = {0.05};
    std::vector<std::filesystem::path> sequence_files = {data("exp_01.fasta")};
    ibf_args.path_out = "Test_";
    std::vector<uint8_t> cutoffs{};
    ibf(sequence_files, ibf_args, minimiser_args, fpr, cutoffs);

    app_test_result result = execute_app("estimate -i ", "Test_", data("mini_gen.fasta"));
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(estimate_options_test, with_argument_normalization_method)
{
    estimate_ibf_arguments ibf_args{};
    minimiser_file_input_arguments minimiser_args{};
    ibf_args.expression_thresholds = {1, 2};
    std::vector<double> fpr = {0.05};
    std::vector<std::filesystem::path> sequence_files = {data("exp_01.fasta")};
    ibf_args.path_out = "Test_";
    std::vector<uint8_t> cutoffs{};
    ibf(sequence_files, ibf_args, minimiser_args, fpr, cutoffs);

    app_test_result result = execute_app("estimate -m -i ", "Test_", data("mini_gen.fasta"));
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(estimate_options_test, with_argument_out)
{
    estimate_ibf_arguments ibf_args{};
    minimiser_file_input_arguments minimiser_args{};
    ibf_args.expression_thresholds = {1, 2};
    std::vector<double> fpr = {0.05};
    std::vector<std::filesystem::path> sequence_files = {data("exp_01.fasta")};
    ibf_args.path_out = "Test_";
    std::vector<uint8_t> cutoffs{};
    ibf(sequence_files, ibf_args, minimiser_args, fpr, cutoffs);

    app_test_result result = execute_app("estimate -o ", "expressions.out", "-i ", "Test_", data("mini_gen.fasta"));
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}
