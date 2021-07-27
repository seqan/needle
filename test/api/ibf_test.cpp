#include <gtest/gtest.h>
#include <iostream>

#include <seqan3/test/expect_range_eq.hpp>

#include "ibf.h"
#include "minimiser.h"

#ifndef DATA_INPUT_DIR
#  define DATA_INPUT_DIR @DATA_INPUT_DIR@
#endif

using seqan3::operator""_shape;
std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory


void initialization_args(arguments & args, ibf_arguments & ibf_args)
{
    args.compressed = true;
    args.k = 4;
    args.shape = seqan3::ungapped{args.k};
    args.w_size = seqan3::window_size{4};
    args.s = seqan3::seed{0};
    args.path_out = tmp_dir/"Test_";
    ibf_args.fpr = {0.05};
}

TEST(ibf, given_expression_levels)
{
    arguments args{};
    ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    initialization_args(args, ibf_args);
    ibf_args.expression_levels = {1, 2};
    std::vector<std::filesystem::path> sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};

    std::vector<uint16_t> expected{1, 2};

    std::vector<uint16_t> medians = ibf(sequence_files, args, ibf_args, minimiser_args);

    EXPECT_EQ(expected, medians);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;

    load_ibf(ibf, tmp_dir/"Test_IBF_1");
    auto agent = ibf.membership_agent();

    std::vector<bool> expected_result(1, 0);
    auto & res = agent.bulk_contains(2);
    EXPECT_RANGE_EQ(expected_result,  res);
    expected_result[0] = 1;
    auto & res2 = agent.bulk_contains(24);
    EXPECT_RANGE_EQ(expected_result,  res2);
    std::filesystem::remove(tmp_dir/"Test_IBF_1");
    std::filesystem::remove(tmp_dir/"Test_IBF_2");
}

TEST(ibf, given_expression_levels_genome_file)
{
    arguments args{};
    ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    initialization_args(args, ibf_args);
    ibf_args.expression_levels = {1, 2};
    minimiser_args.include_file = std::string(DATA_INPUT_DIR) + "mini_example.fasta";
    std::vector<std::filesystem::path> sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};

    std::vector<uint16_t> expected{1, 2};

    std::vector<uint16_t> medians = ibf(sequence_files, args, ibf_args, minimiser_args);

    EXPECT_EQ(expected, medians);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;

    load_ibf(ibf, tmp_dir/"Test_IBF_1");
    auto agent = ibf.membership_agent();

    std::vector<bool> expected_result(1, 0);
    auto & res = agent.bulk_contains(2);
    EXPECT_RANGE_EQ(expected_result,  res);
    expected_result[0] = 1;
    auto & res2 = agent.bulk_contains(24);
    EXPECT_RANGE_EQ(expected_result,  res2);
    std::filesystem::remove(tmp_dir/"Test_IBF_1");
    std::filesystem::remove(tmp_dir/"Test_IBF_2");
}

TEST(ibf, no_given_expression_levels)
{
    arguments args{};
    ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    initialization_args(args, ibf_args);
    ibf_args.number_expression_levels = 2;
    std::vector<std::filesystem::path> sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};

    std::vector<uint16_t> expected{};

    std::vector<uint16_t> medians = ibf(sequence_files, args, ibf_args, minimiser_args);

    EXPECT_EQ(expected, medians);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    load_ibf(ibf, tmp_dir/"Test_IBF_Level_0");
    auto agent = ibf.membership_agent();

    std::vector<bool> expected_result(1, 0);
    auto & res = agent.bulk_contains(2);
    EXPECT_RANGE_EQ(expected_result,  res);
    expected_result[0] = 1;
    auto & res2 = agent.bulk_contains(97);
    EXPECT_RANGE_EQ(expected_result,  res2);
    std::filesystem::remove(tmp_dir/"Test_IBF_Level_0");
}

TEST(ibf, throws)
{
    arguments args{};
    ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    initialization_args(args, ibf_args);
    std::vector<std::filesystem::path> sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};

    EXPECT_THROW(ibf(sequence_files, args, ibf_args, minimiser_args), std::invalid_argument);

    ibf_args.number_expression_levels = 0;
    ibf_args.fpr = {};
    EXPECT_THROW(ibf(sequence_files, args, ibf_args, minimiser_args), std::invalid_argument);

    ibf_args.fpr = {0.05};
    EXPECT_THROW(ibf(sequence_files, args, ibf_args, minimiser_args), std::invalid_argument);
}
