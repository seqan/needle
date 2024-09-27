// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <string> // strings

#include "../app_test.hpp"
#include "ibf.hpp"
#include "shared.hpp"

struct delete_options_test : public app_test
{};

TEST_F(delete_options_test, delete_no_options)
{
    app_test_result result = execute_app("delete");
    std::string expected{"needle-delete - Delete experiments specified by their position from the Needle index.\n"
                         "=====================================================================================\n"
                         "    Try -h or --help for more information.\n"};
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(delete_options_test, with_argument)
{
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    ibf_args.expression_thresholds = {1, 2};
    std::vector<double> fpr = {0.05};
    std::vector<std::filesystem::path> sequence_files = {data("exp_01.fasta"), data("exp_01.fasta")};
    ibf_args.path_out = "Test_";
    std::vector<uint8_t> cutoffs{};
    ibf(sequence_files, ibf_args, minimiser_args, fpr, cutoffs);

    app_test_result result = execute_app("delete -i ", "Test_ ", "0");
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}
