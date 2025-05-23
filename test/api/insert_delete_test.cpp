// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <iostream>

#include <seqan3/test/expect_range_eq.hpp>

#include "../app_test.hpp"
#include "delete.hpp"
#include "ibf.hpp"
#include "insert.hpp"
#include "misc/read_levels.hpp"
#include "shared.hpp"

// To prevent issues when running multiple CLI tests in parallel, give each CLI test unique names:
struct delete_test : public app_test
{
    void initialization_args(estimate_ibf_arguments & args)
    {
        args.k = 4;
        args.shape = seqan3::ungapped{args.k};
        args.w_size = seqan3::window_size{4};
        args.s = seqan3::seed{0};
    }
};

struct insert_test : public delete_test
{};

TEST_F(delete_test, no_given_thresholds)
{
    std::vector<double> fpr = {0.05};
    std::vector<uint8_t> cutoffs_delete{0, 0};
    estimate_ibf_arguments ibf_args_delete{};
    minimiser_file_input_arguments minimiser_args_delete{};
    initialization_args(ibf_args_delete);
    ibf_args_delete.number_expression_thresholds = 2;
    minimiser_args_delete.experiment_names = false;
    ibf_args_delete.path_out = "IBF_delete_Exp_";
    std::vector<std::filesystem::path> sequence_files_delete = {data("mini_example.fasta"), data("mini_example.fasta")};
    ibf(sequence_files_delete, ibf_args_delete, minimiser_args_delete, fpr, cutoffs_delete);
    seqan::hibf::interleaved_bloom_filter ibf{};
    load_ibf(ibf, "IBF_delete_Exp_IBF_Level_0");
    seqan::hibf::interleaved_bloom_filter ibf_0{seqan::hibf::bin_count{2u},
                                                seqan::hibf::bin_size{ibf.bin_size()},
                                                seqan::hibf::hash_function_count{1u}};
    load_ibf(ibf, "IBF_delete_Exp_IBF_Level_1");
    seqan::hibf::interleaved_bloom_filter ibf_1{seqan::hibf::bin_count{2u},
                                                seqan::hibf::bin_size{ibf.bin_size()},
                                                seqan::hibf::hash_function_count{1u}};

    delete_bin({0, 1}, ibf_args_delete, "IBF_delete_Exp_", true);

    seqan::hibf::interleaved_bloom_filter ibf_delete{};

    load_ibf(ibf_delete, "IBF_delete_Exp_IBF_Level_0");
    EXPECT_TRUE((ibf_0 == ibf_delete));

    load_ibf(ibf_delete, "IBF_delete_Exp_IBF_Level_1");
    EXPECT_TRUE((ibf_1 == ibf_delete));
}

TEST_F(insert_test, ibf)
{
    estimate_ibf_arguments ibf_args{};
    minimiser_file_input_arguments minimiser_args{};
    initialization_args(ibf_args);
    ibf_args.path_out = "IBF_True_Exp_";
    ibf_args.expression_thresholds = {1, 2};
    minimiser_args.experiment_names = false;
    std::vector<std::filesystem::path> sequence_files = {data("mini_example.fasta"), data("mini_example.fasta")};
    std::vector<double> fpr = {0.05};
    std::vector<uint8_t> cutoffs{0, 0};

    ibf(sequence_files, ibf_args, minimiser_args, fpr, cutoffs);

    std::vector<uint8_t> cutoffs_insert{0};
    estimate_ibf_arguments ibf_args_insert{};
    minimiser_file_input_arguments minimiser_args_insert{};
    initialization_args(ibf_args_insert);
    ibf_args_insert.path_out = "IBF_True_Exp_";
    ibf_args_insert.expression_thresholds = {1, 2};
    minimiser_args_insert.experiment_names = false;
    ibf_args_insert.path_out = "IBF_Insert_Exp_";
    std::vector<std::filesystem::path> sequence_files_insert = {data("mini_example.fasta")};
    ibf(sequence_files_insert, ibf_args_insert, minimiser_args_insert, fpr, cutoffs_insert);

    insert(sequence_files_insert, ibf_args_insert, minimiser_args_insert, cutoffs_insert, "", "IBF_Insert_Exp_", false);

    seqan::hibf::interleaved_bloom_filter ibf;
    seqan::hibf::interleaved_bloom_filter ibf_insert;

    load_ibf(ibf, "IBF_True_Exp_IBF_1");
    load_ibf(ibf_insert, "IBF_Insert_Exp_IBF_1");
    EXPECT_TRUE((ibf == ibf_insert));

    load_ibf(ibf, "IBF_True_Exp_IBF_2");
    load_ibf(ibf_insert, "IBF_Insert_Exp_IBF_2");
    EXPECT_TRUE((ibf == ibf_insert));

    std::vector<std::vector<double>> fpr_ibf{};
    read_levels<double>(fpr_ibf, "IBF_True_Exp_IBF_FPRs.fprs");
    std::vector<std::vector<double>> fpr_insert{};
    read_levels<double>(fpr_insert, "IBF_Insert_Exp_IBF_FPRs.fprs");
    EXPECT_EQ(fpr_ibf, fpr_insert);
}

TEST_F(insert_test, ibf_no_given_thresholds)
{
    estimate_ibf_arguments ibf_args{};
    minimiser_file_input_arguments minimiser_args{};
    initialization_args(ibf_args);
    ibf_args.path_out = "IBF_True_Exp_";
    ibf_args.number_expression_thresholds = 2;
    minimiser_args.experiment_names = false;
    std::vector<std::filesystem::path> sequence_files = {data("mini_example.fasta"), data("mini_example.fasta")};
    std::vector<double> fpr = {0.05};

    std::vector<uint16_t> expected{};
    std::vector<uint8_t> cutoffs{0, 0};

    std::vector<uint16_t> medians = ibf(sequence_files, ibf_args, minimiser_args, fpr, cutoffs);

    std::vector<uint8_t> cutoffs_insert{0};
    estimate_ibf_arguments ibf_args_insert{};
    minimiser_file_input_arguments minimiser_args_insert{};
    initialization_args(ibf_args_insert);
    ibf_args_insert.number_expression_thresholds = 2;
    minimiser_args_insert.experiment_names = false;
    ibf_args_insert.path_out = "IBF_Insert_Exp_";
    std::vector<std::filesystem::path> sequence_files_insert = {data("mini_example.fasta")};
    std::vector<uint16_t> medians_insert =
        ibf(sequence_files_insert, ibf_args_insert, minimiser_args_insert, fpr, cutoffs_insert);
    insert(sequence_files_insert, ibf_args_insert, minimiser_args_insert, cutoffs_insert, "", "IBF_Insert_Exp_", true);
    seqan::hibf::interleaved_bloom_filter ibf;
    seqan::hibf::interleaved_bloom_filter ibf_insert;

    load_ibf(ibf, "IBF_True_Exp_IBF_Level_0");
    load_ibf(ibf_insert, "IBF_Insert_Exp_IBF_Level_0");
    EXPECT_TRUE((ibf == ibf_insert));

    load_ibf(ibf, "IBF_True_Exp_IBF_Level_1");
    load_ibf(ibf_insert, "IBF_Insert_Exp_IBF_Level_1");
    EXPECT_TRUE((ibf == ibf_insert));

    std::vector<std::vector<uint16_t>> expressions_ibf{};
    read_levels<uint16_t>(expressions_ibf, "IBF_True_Exp_IBF_Levels.levels");
    std::vector<std::vector<uint16_t>> expressions_insert{};
    read_levels<uint16_t>(expressions_insert, "IBF_Insert_Exp_IBF_Levels.levels");
    EXPECT_EQ(expressions_ibf, expressions_insert);

    std::vector<std::vector<double>> fpr_ibf{};
    read_levels<double>(fpr_ibf, "IBF_True_Exp_IBF_FPRs.fprs");
    std::vector<std::vector<double>> fpr_insert{};
    read_levels<double>(fpr_insert, "IBF_Insert_Exp_IBF_FPRs.fprs");
    EXPECT_EQ(fpr_ibf, fpr_insert);
}

TEST_F(insert_test, ibf_delete)
{
    estimate_ibf_arguments ibf_args{};
    minimiser_file_input_arguments minimiser_args{};
    initialization_args(ibf_args);
    ibf_args.path_out = "IBF_True_Exp_";
    ibf_args.expression_thresholds = {1, 2};
    minimiser_args.experiment_names = false;
    std::vector<std::filesystem::path> sequence_files = {data("mini_example.fasta"),
                                                         data("mini_example2.fasta"),
                                                         data("mini_example.fasta")};
    std::vector<double> fpr = {0.05, 0.05};
    std::vector<uint8_t> cutoffs{0, 0, 0};

    ibf(sequence_files, ibf_args, minimiser_args, fpr, cutoffs);

    std::vector<uint8_t> cutoffs_insert{0};
    estimate_ibf_arguments ibf_args_insert{};
    minimiser_file_input_arguments minimiser_args_insert{};
    initialization_args(ibf_args_insert);
    ibf_args_insert.path_out = "IBF_Insert_Exp_";
    ibf_args_insert.expression_thresholds = {1, 2};
    minimiser_args_insert.experiment_names = false;
    std::vector<std::filesystem::path> sequence_files_insert = {data("mini_example2.fasta")};
    std::vector<std::filesystem::path> sequence_files_test = {data("mini_example.fasta"),
                                                              data("mini_example2.fasta"),
                                                              data("mini_example.fasta")};

    ibf(sequence_files_test, ibf_args_insert, minimiser_args, fpr, cutoffs);
    delete_bin({1}, ibf_args_insert, ibf_args_insert.path_out, false);
    insert(sequence_files_insert, ibf_args_insert, minimiser_args_insert, cutoffs_insert, "", "IBF_Insert_Exp_", false);

    seqan::hibf::interleaved_bloom_filter ibf;
    seqan::hibf::interleaved_bloom_filter ibf_insert;

    load_ibf(ibf, "IBF_True_Exp_IBF_1");
    load_ibf(ibf_insert, "IBF_Insert_Exp_IBF_1");
    EXPECT_TRUE((ibf == ibf_insert));

    load_ibf(ibf, "IBF_True_Exp_IBF_2");
    load_ibf(ibf_insert, "IBF_Insert_Exp_IBF_2");
    EXPECT_TRUE((ibf == ibf_insert));

    std::vector<std::vector<double>> fpr_ibf{};
    read_levels<double>(fpr_ibf, "IBF_True_Exp_IBF_FPRs.fprs");
    std::vector<std::vector<double>> fpr_insert{};
    read_levels<double>(fpr_insert, "IBF_Insert_Exp_IBF_FPRs.fprs");
    EXPECT_EQ(fpr_ibf, fpr_insert);
}

TEST_F(insert_test, ibf_delete_no_given_threshold)
{
    estimate_ibf_arguments ibf_args{};
    minimiser_file_input_arguments minimiser_args{};
    initialization_args(ibf_args);
    ibf_args.path_out = "IBF_True_Exp_";
    ibf_args.number_expression_thresholds = 2;
    minimiser_args.experiment_names = false;
    std::vector<std::filesystem::path> sequence_files = {data("mini_example.fasta"),
                                                         data("mini_example2.fasta"),
                                                         data("mini_example.fasta")};
    std::vector<double> fpr = {0.05, 0.05};
    std::vector<uint8_t> cutoffs{0, 0, 0};

    ibf(sequence_files, ibf_args, minimiser_args, fpr, cutoffs);

    std::vector<uint8_t> cutoffs_insert{0};
    estimate_ibf_arguments ibf_args_insert{};
    minimiser_file_input_arguments minimiser_args_insert{};
    initialization_args(ibf_args_insert);
    ibf_args_insert.path_out = "IBF_Insert_Exp_";
    ibf_args_insert.number_expression_thresholds = 2;
    minimiser_args_insert.experiment_names = false;
    std::vector<std::filesystem::path> sequence_files_insert = {data("mini_example2.fasta")};
    std::vector<std::filesystem::path> sequence_files_test = {data("mini_example.fasta"),
                                                              data("mini_example2.fasta"),
                                                              data("mini_example.fasta")};

    ibf(sequence_files_test, ibf_args_insert, minimiser_args, fpr, cutoffs);
    delete_bin({1}, ibf_args_insert, ibf_args_insert.path_out, true);
    insert(sequence_files_insert, ibf_args_insert, minimiser_args_insert, cutoffs_insert, "", "IBF_Insert_Exp_", true);

    seqan::hibf::interleaved_bloom_filter ibf;
    seqan::hibf::interleaved_bloom_filter ibf_insert;

    load_ibf(ibf, "IBF_True_Exp_IBF_Level_0");

    load_ibf(ibf_insert, "IBF_Insert_Exp_IBF_Level_0");
    EXPECT_TRUE((ibf == ibf_insert));

    load_ibf(ibf, "IBF_True_Exp_IBF_Level_1");
    load_ibf(ibf_insert, "IBF_Insert_Exp_IBF_Level_1");
    EXPECT_TRUE((ibf == ibf_insert));

    std::vector<std::vector<uint16_t>> expressions_ibf{};
    read_levels<uint16_t>(expressions_ibf, "IBF_True_Exp_IBF_Levels.levels");
    std::vector<std::vector<uint16_t>> expressions_insert{};
    read_levels<uint16_t>(expressions_insert, "IBF_Insert_Exp_IBF_Levels.levels");
    EXPECT_EQ(expressions_ibf, expressions_insert);

    std::vector<std::vector<double>> fpr_ibf{};
    read_levels<double>(fpr_ibf, "IBF_True_Exp_IBF_FPRs.fprs");
    std::vector<std::vector<double>> fpr_insert{};
    read_levels<double>(fpr_insert, "IBF_Insert_Exp_IBF_FPRs.fprs");
    EXPECT_EQ(fpr_ibf, fpr_insert);
}

TEST_F(insert_test, ibfmin)
{
    estimate_ibf_arguments ibf_args{};
    initialization_args(ibf_args);
    ibf_args.expression_thresholds = {1, 2};
    std::vector<double> fpr = {0.05, 0.05};
    ibf_args.path_out = "IBFMIN_Test_Given_";
    std::vector<std::filesystem::path> minimiser_file = {data("mini_example.minimiser"),
                                                         data("mini_example.minimiser")};
    ibf(minimiser_file, ibf_args, fpr);

    estimate_ibf_arguments ibf_args_insert{};
    initialization_args(ibf_args_insert);
    ibf_args_insert.expression_thresholds = {1, 2};
    ibf_args_insert.path_out = "IBFMIN_Insert_Given_";
    std::vector<std::filesystem::path> minimiser_file_insert = {data("mini_example.minimiser")};
    ibf(minimiser_file_insert, ibf_args_insert, fpr);
    insert(minimiser_file_insert, ibf_args_insert, "", "IBFMIN_Insert_Given_", false);

    seqan::hibf::interleaved_bloom_filter ibf;
    seqan::hibf::interleaved_bloom_filter ibf_insert;

    load_ibf(ibf, "IBFMIN_Test_Given_IBF_1");
    load_ibf(ibf_insert, "IBFMIN_Insert_Given_IBF_1");
    EXPECT_TRUE((ibf == ibf_insert));

    load_ibf(ibf, "IBFMIN_Test_Given_IBF_2");
    load_ibf(ibf_insert, "IBFMIN_Insert_Given_IBF_2");
    EXPECT_TRUE((ibf == ibf_insert));

    std::vector<std::vector<double>> fpr_ibf{};
    read_levels<double>(fpr_ibf, "IBFMIN_Test_Given_IBF_FPRs.fprs");
    std::vector<std::vector<double>> fpr_insert{};
    read_levels<double>(fpr_insert, "IBFMIN_Insert_Given_IBF_FPRs.fprs");
    EXPECT_EQ(fpr_ibf, fpr_insert);
}

TEST_F(insert_test, ibfmin_delete)
{
    estimate_ibf_arguments ibf_args{};
    initialization_args(ibf_args);
    ibf_args.expression_thresholds = {1, 2};
    std::vector<double> fpr = {0.05, 0.05};
    ibf_args.path_out = "IBFMIN_Test_Given_";
    std::vector<std::filesystem::path> minimiser_file = {data("mini_example.minimiser"),
                                                         data("mini_example.minimiser"),
                                                         data("mini_example.minimiser")};
    ibf(minimiser_file, ibf_args, fpr);

    estimate_ibf_arguments ibf_args_insert{};
    initialization_args(ibf_args_insert);
    ibf_args_insert.expression_thresholds = {1, 2};
    ibf_args_insert.path_out = "IBFMIN_Insert_Given_";
    std::vector<std::filesystem::path> minimiser_file_insert = {data("mini_example.minimiser")};
    ibf(minimiser_file, ibf_args_insert, fpr);
    delete_bin({1}, ibf_args_insert, ibf_args_insert.path_out, false);
    insert(minimiser_file_insert, ibf_args_insert, "", "IBFMIN_Insert_Given_", false);

    seqan::hibf::interleaved_bloom_filter ibf;
    seqan::hibf::interleaved_bloom_filter ibf_insert;

    load_ibf(ibf, "IBFMIN_Test_Given_IBF_1");
    load_ibf(ibf_insert, "IBFMIN_Insert_Given_IBF_1");
    EXPECT_TRUE((ibf == ibf_insert));

    load_ibf(ibf, "IBFMIN_Test_Given_IBF_2");
    load_ibf(ibf_insert, "IBFMIN_Insert_Given_IBF_2");
    EXPECT_TRUE((ibf == ibf_insert));

    std::vector<std::vector<double>> fpr_ibf{};
    read_levels<double>(fpr_ibf, "IBFMIN_Test_Given_IBF_FPRs.fprs");
    std::vector<std::vector<double>> fpr_insert{};
    read_levels<double>(fpr_insert, "IBFMIN_Insert_Given_IBF_FPRs.fprs");
    EXPECT_EQ(fpr_ibf, fpr_insert);
}

TEST_F(insert_test, ibfmin_no_given_thresholds)
{
    estimate_ibf_arguments ibf_args{};
    initialization_args(ibf_args);
    ibf_args.number_expression_thresholds = 2;
    std::vector<double> fpr = {0.05, 0.05};
    ibf_args.path_out = "IBFMIN_Test_Given_";
    std::vector<std::filesystem::path> minimiser_file = {data("mini_example.minimiser"),
                                                         data("mini_example.minimiser")};

    ibf(minimiser_file, ibf_args, fpr);

    estimate_ibf_arguments ibf_args_insert{};
    initialization_args(ibf_args_insert);
    ibf_args_insert.number_expression_thresholds = 2;
    ibf_args_insert.path_out = "IBFMIN_Insert_Given_";
    fpr = {0.05};
    std::vector<std::filesystem::path> minimiser_file_insert = {data("mini_example.minimiser")};
    ibf(minimiser_file_insert, ibf_args_insert, fpr);
    insert(minimiser_file_insert, ibf_args_insert, "", "IBFMIN_Insert_Given_", true);

    seqan::hibf::interleaved_bloom_filter ibf;
    seqan::hibf::interleaved_bloom_filter ibf_insert;

    load_ibf(ibf, "IBFMIN_Test_Given_IBF_Level_0");
    load_ibf(ibf_insert, "IBFMIN_Insert_Given_IBF_Level_0");
    EXPECT_TRUE((ibf == ibf_insert));

    load_ibf(ibf, "IBFMIN_Test_Given_IBF_Level_1");
    load_ibf(ibf_insert, "IBFMIN_Insert_Given_IBF_Level_1");
    EXPECT_TRUE((ibf == ibf_insert));

    std::vector<std::vector<uint16_t>> expressions_ibf{};
    read_levels<uint16_t>(expressions_ibf, "IBFMIN_Test_Given_IBF_Levels.levels");
    std::vector<std::vector<uint16_t>> expressions_insert{};
    read_levels<uint16_t>(expressions_insert, "IBFMIN_Insert_Given_IBF_Levels.levels");
    EXPECT_EQ(expressions_ibf, expressions_insert);

    std::vector<std::vector<double>> fpr_ibf{};
    read_levels<double>(fpr_ibf, "IBFMIN_Test_Given_IBF_FPRs.fprs");
    std::vector<std::vector<double>> fpr_insert{};
    read_levels<double>(fpr_insert, "IBFMIN_Insert_Given_IBF_FPRs.fprs");
    EXPECT_EQ(fpr_ibf, fpr_insert);
}

TEST_F(insert_test, delete_ibfmin_no_given_thresholds)
{
    estimate_ibf_arguments ibf_args{};
    initialization_args(ibf_args);
    ibf_args.number_expression_thresholds = 2;
    std::vector<double> fpr = {0.05, 0.05};
    ibf_args.path_out = "IBFMIN_Test_Given_Del_";
    std::vector<std::filesystem::path> minimiser_file = {data("mini_example.minimiser"),
                                                         data("mini_example.minimiser"),
                                                         data("mini_example.minimiser")};

    ibf(minimiser_file, ibf_args, fpr);

    estimate_ibf_arguments ibf_args_insert{};
    initialization_args(ibf_args_insert);
    ibf_args_insert.number_expression_thresholds = 2;
    ibf_args_insert.path_out = "IBFMIN_Insert_Given_Del_";
    fpr = {0.05};
    std::vector<std::filesystem::path> minimiser_file_insert = {data("mini_example.minimiser")};
    ibf(minimiser_file, ibf_args_insert, fpr);
    delete_bin({1}, ibf_args_insert, ibf_args_insert.path_out, true);
    insert(minimiser_file_insert, ibf_args_insert, "", "IBFMIN_Insert_Given_Del_", true);

    seqan::hibf::interleaved_bloom_filter ibf;
    seqan::hibf::interleaved_bloom_filter ibf_insert;

    load_ibf(ibf, "IBFMIN_Test_Given_Del_IBF_Level_0");
    load_ibf(ibf_insert, "IBFMIN_Insert_Given_Del_IBF_Level_0");
    EXPECT_TRUE((ibf == ibf_insert));

    load_ibf(ibf, "IBFMIN_Test_Given_Del_IBF_Level_1");
    load_ibf(ibf_insert, "IBFMIN_Insert_Given_Del_IBF_Level_1");
    EXPECT_TRUE((ibf == ibf_insert));

    std::vector<std::vector<uint16_t>> expressions_ibf{};
    read_levels<uint16_t>(expressions_ibf, "IBFMIN_Test_Given_Del_IBF_Levels.levels");
    std::vector<std::vector<uint16_t>> expressions_insert{};
    read_levels<uint16_t>(expressions_insert, "IBFMIN_Insert_Given_Del_IBF_Levels.levels");
    EXPECT_EQ(expressions_ibf, expressions_insert);

    std::vector<std::vector<double>> fpr_ibf{};
    read_levels<double>(fpr_ibf, "IBFMIN_Test_Given_Del_IBF_FPRs.fprs");
    std::vector<std::vector<double>> fpr_insert{};
    read_levels<double>(fpr_insert, "IBFMIN_Insert_Given_Del_IBF_FPRs.fprs");
    EXPECT_EQ(fpr_ibf, fpr_insert);
}
