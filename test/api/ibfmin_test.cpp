#include <gtest/gtest.h>
#include <iostream>

#include <seqan3/test/expect_range_eq.hpp>

#include "ibf.h"
#include "shared.h"

#ifndef DATA_INPUT_DIR
#  define DATA_INPUT_DIR @DATA_INPUT_DIR@
#endif

using seqan3::operator""_shape;
std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory

void initialization_args(estimate_ibf_arguments & args)
{
    args.compressed = true;
    args.k = 4;
    args.shape = seqan3::ungapped{args.k};
    args.w_size = seqan3::window_size{4};
    args.s = seqan3::seed{0};
    args.path_out = tmp_dir/"IBFMIN_Test_";
}

TEST(ibfmin, given_expression_levels)
{
    estimate_ibf_arguments ibf_args{};
    initialization_args(ibf_args);
    ibf_args.expression_levels = {1, 2};
    std::vector<double> fpr = {0.05, 0.05};
    ibf_args.path_out = tmp_dir/"IBFMIN_Test_Given_";
    std::vector<std::filesystem::path> minimiser_file = {std::string(DATA_INPUT_DIR) + "mini_example.minimiser"};

    std::vector<uint16_t> expected{1, 2};

    std::vector<uint16_t> medians = ibf(minimiser_file, ibf_args, fpr);

    EXPECT_EQ(expected, medians);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    if (std::filesystem::exists(tmp_dir/"IBFMIN_Test_Given_IBF_1"))
    {
        load_ibf(ibf, tmp_dir/"IBFMIN_Test_Given_IBF_1");
        auto agent = ibf.membership_agent();

        std::vector<bool> expected_result(1, 0);
        auto & res = agent.bulk_contains(97);
        EXPECT_RANGE_EQ(expected_result,  res);
        expected_result[0] = 1;
        auto & res2 = agent.bulk_contains(24);
        EXPECT_RANGE_EQ(expected_result,  res2);
    }

    std::filesystem::remove(tmp_dir/"IBFMIN_Test_Given_IBF_1");
    std::filesystem::remove(tmp_dir/"IBFMIN_Test_Given_IBF_2");
}

#if defined(__GNUC__) && ((__GNUC___ == 10 && __cplusplus == 201703L) || (__GNUC__ <10))
TEST(ibfmin, given_expression_levels_multiple_threads)
{
    estimate_ibf_arguments ibf_args{};
    initialization_args(ibf_args);
    ibf_args.expression_levels = {1, 2};
    std::vector<double> fpr = {0.05, 0.05};
    ibf_args.threads = 2;
    ibf_args.path_out = tmp_dir/"IBFMIN_Test_Multiple_";
    std::vector<std::filesystem::path> minimiser_file{};
    minimiser_file.assign(16, std::string(DATA_INPUT_DIR) + "mini_example.minimiser");

    std::vector<uint16_t> expected{1, 2};

    std::vector<uint16_t> medians = ibf(minimiser_file, ibf_args, fpr);

    EXPECT_EQ(expected, medians);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    load_ibf(ibf, tmp_dir/"IBFMIN_Test_Multiple_IBF_1");
    auto agent = ibf.membership_agent();

    std::vector<bool> expected_result(16, 0);
    auto & res = agent.bulk_contains(97);
    EXPECT_RANGE_EQ(expected_result,  res);
    std::vector<bool> expected_result2(16, 1);
    auto & res2 = agent.bulk_contains(24);
    EXPECT_RANGE_EQ(expected_result2,  res2);
    std::filesystem::remove(tmp_dir/"IBFMIN_Test_Multiple_IBF_1");
    std::filesystem::remove(tmp_dir/"IBFMIN_Test_Multiple_IBF_2");
}
#endif

TEST(ibfmin, no_given_expression_levels)
{
    estimate_ibf_arguments ibf_args{};
    initialization_args(ibf_args);
    ibf_args.number_expression_levels = 2;
    std::vector<double> fpr = {0.0025, 0.0025};
    std::vector<std::filesystem::path> minimiser_file = {std::string(DATA_INPUT_DIR) + "mini_example.minimiser"};

    std::vector<uint16_t> expected{};

    std::vector<uint16_t> medians = ibf(minimiser_file, ibf_args, fpr);

    EXPECT_EQ(expected, medians);

    if (std::filesystem::exists(tmp_dir/"IBFMIN_Test_IBF_Level_0"))
    {
        seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
        load_ibf(ibf, tmp_dir/"IBFMIN_Test_IBF_Level_0");
        auto agent = ibf.membership_agent();

        std::vector<bool> expected_result(1, 0);
        auto & res = agent.bulk_contains(2);
        EXPECT_RANGE_EQ(expected_result,  res);
        expected_result[0] = 1;
        auto & res2 = agent.bulk_contains(97);
        EXPECT_RANGE_EQ(expected_result,  res2);
    }

    std::filesystem::remove(tmp_dir/"IBFMIN_Test_IBF_Level_0");
    std::filesystem::remove(tmp_dir/"IBFMIN_Test_IBF_Level_1");
    std::filesystem::remove(tmp_dir/"IBFMIN_Test_IBF_Levels.levels");
}

TEST(ibfmin, expression_levels_by_genome)
{
    estimate_ibf_arguments ibf_args{};
    initialization_args(ibf_args);
    ibf_args.number_expression_levels = 1;
    std::vector<double> fpr = {0.05};
    std::vector<std::filesystem::path> minimiser_file = {std::string(DATA_INPUT_DIR) + "mini_example.minimiser"};

    std::vector<uint16_t> expected{};

    std::vector<uint16_t> medians = ibf(minimiser_file, ibf_args, fpr, std::string(DATA_INPUT_DIR) + "mini_gen.fasta");

    EXPECT_EQ(expected, medians);

    if (std::filesystem::exists(tmp_dir/"IBFMIN_Test_IBF_Level_0"))
    {
        seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
        load_ibf(ibf, tmp_dir/"IBFMIN_Test_IBF_Level_0");
        auto agent = ibf.membership_agent();

        std::vector<bool> expected_result(1, 0);
        auto & res = agent.bulk_contains(0);
        EXPECT_RANGE_EQ(expected_result,  res);
        expected_result[0] = 1;
        auto & res2 = agent.bulk_contains(97);
        EXPECT_RANGE_EQ(expected_result,  res2);
    }

    std::filesystem::remove(tmp_dir/"IBFMIN_Test_IBF_Level_0");
    std::filesystem::remove(tmp_dir/"IBFMIN_Test_IBF_Levels.levels");
}

#if defined(__GNUC__) && ((__GNUC___ == 10 && __cplusplus == 201703L) || (__GNUC__ <10))
TEST(ibfmin, no_given_expression_levels_multiple_threads)
{
    estimate_ibf_arguments ibf_args{};
    initialization_args(ibf_args);
    ibf_args.number_expression_levels = 2;
    std::vector<double> fpr = {0.0025, 0.0025};
    ibf_args.threads = 2;
    ibf_args.path_out = tmp_dir/"IBFMIN_Test_Multiple_";
    std::vector<std::filesystem::path> minimiser_file{};
    minimiser_file.assign(128, std::string(DATA_INPUT_DIR) + "mini_example.minimiser");

    std::vector<uint16_t> expected{};

    std::vector<uint16_t> medians = ibf(minimiser_file, ibf_args, fpr);

    EXPECT_EQ(expected, medians);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    load_ibf(ibf, tmp_dir/"IBFMIN_Test_Multiple_IBF_Level_0");
    auto agent = ibf.membership_agent();

    std::vector<bool> expected_result(128, 0);
    auto & res = agent.bulk_contains(2);
    EXPECT_RANGE_EQ(expected_result,  res);
    std::vector<bool> expected_result2(128, 1);
    auto & res2 = agent.bulk_contains(97);
    EXPECT_RANGE_EQ(expected_result2,  res2);
    std::filesystem::remove(tmp_dir/"IBFMIN_Test_Multiple_IBF_Level_0");
    std::filesystem::remove(tmp_dir/"IBFMIN_Test_Multiple_IBF_Level_1");
    std::filesystem::remove(tmp_dir/"IBFMIN_Test_Multiple_IBF_Levels.levels");
}
#endif
