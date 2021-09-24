#include <gtest/gtest.h>
#include <iostream>

#include <seqan3/test/expect_range_eq.hpp>

#include "ibf.h"
#include "shared.h"

#ifndef DATA_INPUT_DIR
#  define DATA_INPUT_DIR @DATA_INPUT_DIR@
#endif

using seqan3::operator""_shape;


void initialization_args(estimate_ibf_arguments & args)
{
    args.compressed = true;
    args.k = 4;
    args.shape = seqan3::ungapped{args.k};
    args.w_size = seqan3::window_size{4};
    args.s = seqan3::seed{0};
}

TEST(ibf, given_expression_thresholds)
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    initialization_args(ibf_args);
    ibf_args.path_out = tmp_dir/"IBF_Test_Exp_";
    ibf_args.expression_thresholds = {1, 2};
    minimiser_args.experiment_names = true;
    std::vector<std::filesystem::path> sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};
    std::vector<double> fpr = {0.05};

    std::vector<uint16_t> expected{1, 2};

    std::vector<uint16_t> medians = ibf(sequence_files, ibf_args, minimiser_args, fpr);

    EXPECT_EQ(expected, medians);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    if (std::filesystem::exists(tmp_dir/"IBF_Test_Exp_IBF_1"))
    {
        load_ibf(ibf, tmp_dir/"IBF_Test_Exp_IBF_1");
        auto agent = ibf.membership_agent();

        std::vector<bool> expected_result(1, 0);
        auto & res = agent.bulk_contains(2);
        EXPECT_RANGE_EQ(expected_result,  res);
        expected_result[0] = 1;
        auto & res2 = agent.bulk_contains(24);
        EXPECT_RANGE_EQ(expected_result,  res2);
    }
    std::filesystem::remove(tmp_dir/"IBF_Test_Exp_IBF_1");
    std::filesystem::remove(tmp_dir/"IBF_Test_Exp_IBF_2");

    estimate_ibf_arguments args{};
    if (std::filesystem::exists(tmp_dir/"IBF_Test_Exp_IBF_Data"))
    {
        load_args(args, tmp_dir/"IBF_Test_Exp_IBF_Data");
        EXPECT_EQ(4, args.k);
        EXPECT_EQ(4, args.w_size.get());
        EXPECT_EQ(0, args.s.get());
        EXPECT_EQ(15, args.shape.to_ulong());
        EXPECT_EQ(true, args.compressed);
        EXPECT_EQ(ibf_args.number_expression_thresholds, args.number_expression_thresholds);
        EXPECT_RANGE_EQ(ibf_args.expression_thresholds, args.expression_thresholds);
        EXPECT_EQ(ibf_args.samplewise, args.samplewise);
    }
    std::filesystem::remove(tmp_dir/"IBF_Test_Exp_IBF_Data");
    std::filesystem::remove(tmp_dir/"IBF_Test_Exp_Test_Stored_Files.txt");
}

TEST(ibf, given_expression_thresholds_include_file)
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    initialization_args(ibf_args);
    ibf_args.path_out = tmp_dir/"IBF_Test_Include_";
    ibf_args.expression_thresholds = {1, 2};
    minimiser_args.include_file = std::string(DATA_INPUT_DIR) + "mini_example.fasta";
    std::vector<std::filesystem::path> sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};
    std::vector<double> fpr = {0.05};

    std::vector<uint16_t> expected{1, 2};

    std::vector<uint16_t> medians = ibf(sequence_files, ibf_args, minimiser_args, fpr);

    EXPECT_EQ(expected, medians);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    if (std::filesystem::exists(tmp_dir/"IBF_Test_Include_IBF_1"))
    {
        load_ibf(ibf, tmp_dir/"IBF_Test_Include_IBF_1");
        auto agent = ibf.membership_agent();

        std::vector<bool> expected_result(1, 0);
        auto & res = agent.bulk_contains(2);
        EXPECT_RANGE_EQ(expected_result,  res);
        expected_result[0] = 1;
        auto & res2 = agent.bulk_contains(24);
        EXPECT_RANGE_EQ(expected_result,  res2);
    }
    std::filesystem::remove(tmp_dir/"IBF_Test_Include_IBF_1");
    std::filesystem::remove(tmp_dir/"IBF_Test_Include_IBF_2");
    std::filesystem::remove(tmp_dir/"IBF_Test_Include_IBF_Data");
}

TEST(ibf, given_expression_thresholds_exclude_file)
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    initialization_args(ibf_args);
    ibf_args.path_out = tmp_dir/"IBF_Test_Exclude_";
    ibf_args.expression_thresholds = {1, 2};
    minimiser_args.exclude_file = std::string(DATA_INPUT_DIR) + "mini_gen.fasta";
    std::vector<std::filesystem::path> sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};
    std::vector<double> fpr = {0.05};

    std::vector<uint16_t> expected{1, 2};

    std::vector<uint16_t> medians = ibf(sequence_files, ibf_args, minimiser_args, fpr);

    EXPECT_EQ(expected, medians);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    if (std::filesystem::exists(tmp_dir/"IBF_Test_Exclude_IBF_1"))
    {
        load_ibf(ibf, tmp_dir/"IBF_Test_Exclude_IBF_1");
        auto agent = ibf.membership_agent();

        std::vector<bool> expected_result(1, 0);
        auto & res = agent.bulk_contains(2);
        EXPECT_RANGE_EQ(expected_result,  res);
        expected_result[0] = 1;
        auto & res2 = agent.bulk_contains(24);
        EXPECT_RANGE_EQ(expected_result,  res2);
    }
    std::filesystem::remove(tmp_dir/"IBF_Test_Exclude_IBF_1");
    std::filesystem::remove(tmp_dir/"IBF_Test_Exclude_IBF_2");
    std::filesystem::remove(tmp_dir/"IBF_Test_Exclude_IBF_Data");
}

TEST(ibf, no_given_expression_thresholds)
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    initialization_args(ibf_args);
    ibf_args.path_out = tmp_dir/"IBF_Test_";
    ibf_args.number_expression_thresholds = 2;
    std::vector<std::filesystem::path> sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};
    std::vector<double> fpr = {0.05};

    std::vector<uint16_t> expected{};

    std::vector<uint16_t> medians = ibf(sequence_files, ibf_args, minimiser_args, fpr);

    EXPECT_EQ(expected, medians);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;

    if (std::filesystem::exists(tmp_dir/"IBF_Test_IBF_Level_0"))
    {
        load_ibf(ibf, tmp_dir/"IBF_Test_IBF_Level_0");
        auto agent = ibf.membership_agent();

        std::vector<bool> expected_result(1, 0);
        auto & res = agent.bulk_contains(2);
        EXPECT_RANGE_EQ(expected_result,  res);
        expected_result[0] = 1;
        auto & res2 = agent.bulk_contains(97);
        EXPECT_RANGE_EQ(expected_result,  res2);
    }
    std::filesystem::remove(tmp_dir/"IBF_Test_IBF_Level_0");
    std::filesystem::remove(tmp_dir/"IBF_Test_IBF_Level_1");
    std::filesystem::remove(tmp_dir/"IBF_Test_IBF_Levels.levels");
    std::filesystem::remove(tmp_dir/"IBF_Test_IBF_Data");
}

TEST(ibf, expression_thresholds_by_genome)
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    initialization_args(ibf_args);
    ibf_args.path_out = tmp_dir/"IBF_Test_";
    ibf_args.number_expression_thresholds = 1;
    std::vector<std::filesystem::path> sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};
    std::vector<double> fpr = {0.05};

    std::vector<uint16_t> expected{};

    std::vector<uint16_t> medians = ibf(sequence_files, ibf_args, minimiser_args, fpr,
                                        std::string(DATA_INPUT_DIR) + "mini_gen.fasta");

    EXPECT_EQ(expected, medians);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;

    if (std::filesystem::exists(tmp_dir/"IBF_Test_IBF_Level_0"))
    {
        load_ibf(ibf, tmp_dir/"IBF_Test_IBF_Level_0");
        auto agent = ibf.membership_agent();

        std::vector<bool> expected_result(1, 0);
        auto & res = agent.bulk_contains(0);
        EXPECT_RANGE_EQ(expected_result,  res);
        expected_result[0] = 1;
        auto & res2 = agent.bulk_contains(97);
        EXPECT_RANGE_EQ(expected_result,  res2);
    }
    std::filesystem::remove(tmp_dir/"IBF_Test_IBF_Level_0");
    std::filesystem::remove(tmp_dir/"IBF_Test_IBF_Levels.levels");
    std::filesystem::remove(tmp_dir/"IBF_Test_IBF_Data");
}

TEST(ibf, throws)
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    initialization_args(ibf_args);
    ibf_args.path_out = tmp_dir/"IBF_Test_";
    std::vector<std::filesystem::path> sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};
    std::vector<double> fpr = {0.05};

    EXPECT_THROW(ibf(sequence_files, ibf_args, minimiser_args, fpr), std::invalid_argument);

    ibf_args.number_expression_thresholds = 0;
    fpr = {};
    EXPECT_THROW(ibf(sequence_files, ibf_args, minimiser_args, fpr), std::invalid_argument);

    fpr = {0.05};
    EXPECT_THROW(ibf(sequence_files, ibf_args, minimiser_args, fpr), std::invalid_argument);

    ibf_args.expression_thresholds = {10, 1000};
    EXPECT_THROW(ibf(sequence_files, ibf_args, minimiser_args, fpr), std::invalid_argument);
}

TEST(ibf, given_cutoffs)
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    initialization_args(ibf_args);
    ibf_args.path_out = tmp_dir/"IBF_Test_Cut_";
    ibf_args.expression_thresholds = {1, 2};
    minimiser_args.cutoffs = {0};
    std::vector<std::filesystem::path> sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};
    std::vector<double> fpr = {0.05};

    std::vector<uint16_t> expected{1, 2};

    std::vector<uint16_t> medians = ibf(sequence_files, ibf_args, minimiser_args, fpr);

    EXPECT_EQ(expected, medians);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    if (std::filesystem::exists(tmp_dir/"IBF_Test_Cut_IBF_1"))
    {
        load_ibf(ibf, tmp_dir/"IBF_Test_Cut_IBF_1");
        auto agent = ibf.membership_agent();

        std::vector<bool> expected_result(1, 0);
        auto & res = agent.bulk_contains(2);
        EXPECT_RANGE_EQ(expected_result,  res);
        expected_result[0] = 1;
        auto & res2 = agent.bulk_contains(24);
        EXPECT_RANGE_EQ(expected_result,  res2);
    }
    std::filesystem::remove(tmp_dir/"IBF_Test_Cut_IBF_1");
    std::filesystem::remove(tmp_dir/"IBF_Test_Cut_IBF_2");

    estimate_ibf_arguments args{};
    if (std::filesystem::exists(tmp_dir/"IBF_Test_Cut_IBF_Data"))
    {
        load_args(args, tmp_dir/"IBF_Test_Cut_IBF_Data");
        EXPECT_EQ(4, args.k);
        EXPECT_EQ(4, args.w_size.get());
        EXPECT_EQ(0, args.s.get());
        EXPECT_EQ(15, args.shape.to_ulong());
        EXPECT_EQ(true, args.compressed);
        EXPECT_EQ(ibf_args.number_expression_thresholds, args.number_expression_thresholds);
        EXPECT_RANGE_EQ(ibf_args.expression_thresholds, args.expression_thresholds);
        EXPECT_EQ(ibf_args.samplewise, args.samplewise);
    }
    std::filesystem::remove(tmp_dir/"IBF_Test_Cut_IBF_Data");
}
