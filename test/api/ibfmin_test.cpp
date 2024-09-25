#include <gtest/gtest.h>
#include <iostream>

#include <seqan3/test/expect_range_eq.hpp>

#include "ibf.h"
#include "shared.h"
#include "../app_test.hpp"

// To prevent issues when running multiple CLI tests in parallel, give each CLI test unique names:
struct ibfmin_test : public app_test
{
    void initialization_args(estimate_ibf_arguments & args)
    {
        args.compressed = true;
        args.k = 4;
        args.shape = seqan3::ungapped{args.k};
        args.w_size = seqan3::window_size{4};
        args.s = seqan3::seed{0};
        args.path_out = "IBFMIN_Test_";
    }
};

TEST_F(ibfmin_test, given_expression_thresholds)
{
    estimate_ibf_arguments ibf_args{};
    initialization_args(ibf_args);
    ibf_args.expression_thresholds = {1, 2};
    std::vector<double> fpr = {0.05, 0.05};
    ibf_args.path_out = "IBFMIN_Test_Given_";
    std::vector<std::filesystem::path> minimiser_file = {data("mini_example.minimiser")};

    std::vector<uint16_t> expected{1, 2};

    std::vector<uint16_t> medians = ibf(minimiser_file, ibf_args, fpr);

    EXPECT_EQ(expected, medians);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    ASSERT_TRUE(std::filesystem::exists("IBFMIN_Test_Given_IBF_1"));
    {
        load_ibf(ibf, "IBFMIN_Test_Given_IBF_1");
        auto agent = ibf.membership_agent();

        std::vector<bool> expected_result(1, 0);
        auto & res = agent.bulk_contains(97);
        EXPECT_RANGE_EQ(expected_result,  res);
        expected_result[0] = 1;
        auto & res2 = agent.bulk_contains(24);
        EXPECT_RANGE_EQ(expected_result,  res2);
    }

}

TEST_F(ibfmin_test, given_expression_thresholds_multiple_threads)
{
    estimate_ibf_arguments ibf_args{};
    initialization_args(ibf_args);
    ibf_args.expression_thresholds = {1, 2};
    std::vector<double> fpr = {0.05, 0.05};
    ibf_args.threads = 2;
    ibf_args.path_out = "IBFMIN_Test_Multiple_";
    std::vector<std::filesystem::path> minimiser_file{};
    minimiser_file.assign(16, data("mini_example.minimiser"));

    std::vector<uint16_t> expected{1, 2};

    std::vector<uint16_t> medians = ibf(minimiser_file, ibf_args, fpr);

    EXPECT_EQ(expected, medians);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    load_ibf(ibf, "IBFMIN_Test_Multiple_IBF_1");
    auto agent = ibf.membership_agent();

    std::vector<bool> expected_result(16, 0);
    auto & res = agent.bulk_contains(97);
    EXPECT_RANGE_EQ(expected_result,  res);
    std::vector<bool> expected_result2(16, 1);
    auto & res2 = agent.bulk_contains(24);
    EXPECT_RANGE_EQ(expected_result2,  res2);
}

TEST_F(ibfmin_test, no_given_expression_thresholds)
{
    estimate_ibf_arguments ibf_args{};
    initialization_args(ibf_args);
    ibf_args.number_expression_thresholds = 2;
    std::vector<double> fpr = {0.0025, 0.0025};
    std::vector<std::filesystem::path> minimiser_file = {data("mini_example.minimiser")};

    std::vector<uint16_t> expected{};

    std::vector<uint16_t> medians = ibf(minimiser_file, ibf_args, fpr);

    EXPECT_EQ(expected, medians);

    ASSERT_TRUE(std::filesystem::exists("IBFMIN_Test_IBF_Level_0"));
    {
        seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
        load_ibf(ibf, "IBFMIN_Test_IBF_Level_0");
        auto agent = ibf.membership_agent();

        std::vector<bool> expected_result(1, 0);
        auto & res = agent.bulk_contains(2);
        EXPECT_RANGE_EQ(expected_result,  res);
        expected_result[0] = 1;
        auto & res2 = agent.bulk_contains(24);
        EXPECT_RANGE_EQ(expected_result,  res2);
    }

}

TEST_F(ibfmin_test, expression_thresholds_by_genome)
{
    estimate_ibf_arguments ibf_args{};
    initialization_args(ibf_args);
    ibf_args.number_expression_thresholds = 1;
    std::vector<double> fpr = {0.05};
    std::vector<std::filesystem::path> minimiser_file = {data("mini_example.minimiser")};

    std::vector<uint16_t> expected{};

    std::vector<uint16_t> medians = ibf(minimiser_file, ibf_args, fpr, data("mini_gen.fasta"));

    EXPECT_EQ(expected, medians);

    ASSERT_TRUE(std::filesystem::exists("IBFMIN_Test_IBF_Level_0"));
    {
        seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
        load_ibf(ibf, "IBFMIN_Test_IBF_Level_0");
        auto agent = ibf.membership_agent();

        std::vector<bool> expected_result(1, 0);
        auto & res = agent.bulk_contains(2);
        EXPECT_RANGE_EQ(expected_result,  res);
        expected_result[0] = 1;
        auto & res2 = agent.bulk_contains(24);
        EXPECT_RANGE_EQ(expected_result,  res2);
    }

}

TEST_F(ibfmin_test, no_given_expression_thresholds_multiple_threads)
{
    estimate_ibf_arguments ibf_args{};
    initialization_args(ibf_args);
    ibf_args.number_expression_thresholds = 2;
    std::vector<double> fpr = {0.0025, 0.0025};
    ibf_args.threads = 2;
    ibf_args.path_out = "IBFMIN_Test_Multiple_";
    std::vector<std::filesystem::path> minimiser_file{};
    minimiser_file.assign(128, data("mini_example.minimiser"));

    std::vector<uint16_t> expected{};

    std::vector<uint16_t> medians = ibf(minimiser_file, ibf_args, fpr);

    EXPECT_EQ(expected, medians);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    load_ibf(ibf, "IBFMIN_Test_Multiple_IBF_Level_0");
    auto agent = ibf.membership_agent();

    std::vector<bool> expected_result(128, 0);
    auto & res = agent.bulk_contains(2);
    EXPECT_RANGE_EQ(expected_result,  res);
    std::vector<bool> expected_result2(128, 1);
    auto & res2 = agent.bulk_contains(24);
    EXPECT_RANGE_EQ(expected_result2,  res2);
}

TEST_F(ibfmin_test, different_shape)
{
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    initialization_args(ibf_args);
    std::vector<uint8_t> cutoffs = {0};
    ibf_args.shape = seqan3::bin_literal{11};
    ibf_args.expression_thresholds = {1, 2};
    std::vector<double> fpr = {0.05, 0.05};
    ibf_args.path_out = "IBFMIN_Test_Shape_";
    std::vector<std::filesystem::path> sequence_files = {data("mini_example.fasta")};
    minimiser(sequence_files, ibf_args, minimiser_args, cutoffs);
    std::vector<std::filesystem::path> minimiser_file = {"IBFMIN_Test_Shape_mini_example.minimiser"};

    std::vector<uint16_t> expected{1, 2};

    std::vector<uint16_t> medians = ibf(minimiser_file, ibf_args, fpr);

    EXPECT_EQ(expected, medians);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    ASSERT_TRUE(std::filesystem::exists("IBFMIN_Test_Shape_IBF_1"));
    {
        load_ibf(ibf, "IBFMIN_Test_Shape_IBF_1");
        auto agent = ibf.membership_agent();

        std::vector<bool> expected_result(1, 0);
        auto & res = agent.bulk_contains(97);
        EXPECT_RANGE_EQ(expected_result,  res);
        expected_result[0] = 1;
        auto & res2 = agent.bulk_contains(4);
        EXPECT_RANGE_EQ(expected_result,  res2);
    }

}
