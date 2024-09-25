#include <gtest/gtest.h>
#include <iostream>

#include <seqan3/test/expect_range_eq.hpp>

#include "ibf.h"
#include "shared.h"
#include "../app_test.hpp"

// To prevent issues when running multiple CLI tests in parallel, give each CLI test unique names:
struct ibf_test : public app_test
{
    void initialization_args(estimate_ibf_arguments & args)
    {
        args.compressed = true;
        args.k = 4;
        args.shape = seqan3::ungapped{args.k};
        args.w_size = seqan3::window_size{4};
        args.s = seqan3::seed{0};
    }
};

TEST_F(ibf_test, given_expression_thresholds)
{
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    initialization_args(ibf_args);
    ibf_args.path_out = "IBF_Test_Exp_";
    ibf_args.expression_thresholds = {1, 2};
    minimiser_args.experiment_names = true;
    std::vector<std::filesystem::path> sequence_files = {data("mini_example.fasta")};
    std::vector<double> fpr = {0.05};

    std::vector<uint16_t> expected{1, 2};
    std::vector<uint8_t> cutoffs{0};

    std::vector<uint16_t> medians = ibf(sequence_files, ibf_args, minimiser_args, fpr, cutoffs);

    EXPECT_EQ(expected, medians);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    ASSERT_TRUE(std::filesystem::exists("IBF_Test_Exp_IBF_1"));
    {
        load_ibf(ibf, "IBF_Test_Exp_IBF_1");
        auto agent = ibf.membership_agent();

        std::vector<bool> expected_result(1, 0);
        auto & res = agent.bulk_contains(2);
        EXPECT_RANGE_EQ(expected_result,  res);
        expected_result[0] = 1;
        auto & res2 = agent.bulk_contains(24);
        EXPECT_RANGE_EQ(expected_result,  res2);
    }

    ASSERT_TRUE(std::filesystem::exists("IBF_Test_Exp_IBF_Data"));
    {
        estimate_ibf_arguments args{};
        load_args(args, "IBF_Test_Exp_IBF_Data");
        EXPECT_EQ(4, args.k);
        EXPECT_EQ(4, args.w_size.get());
        EXPECT_EQ(0, args.s.get());
        EXPECT_EQ(15, args.shape.to_ulong());
        EXPECT_EQ(true, args.compressed);
        EXPECT_EQ(ibf_args.number_expression_thresholds, args.number_expression_thresholds);
        EXPECT_RANGE_EQ(ibf_args.expression_thresholds, args.expression_thresholds);
        EXPECT_EQ(ibf_args.samplewise, args.samplewise);
    }
}

TEST_F(ibf_test, given_expression_thresholds_include_file)
{
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    initialization_args(ibf_args);
    ibf_args.path_out = "IBF_Test_Include_";
    ibf_args.expression_thresholds = {1, 2};
    minimiser_args.include_file = data("mini_example.fasta");
    std::vector<std::filesystem::path> sequence_files = {data("mini_example.fasta")};
    std::vector<double> fpr = {0.05};

    std::vector<uint16_t> expected{1, 2};
    std::vector<uint8_t> cutoffs{0};

    std::vector<uint16_t> medians = ibf(sequence_files, ibf_args, minimiser_args, fpr, cutoffs);

    EXPECT_EQ(expected, medians);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    ASSERT_TRUE(std::filesystem::exists("IBF_Test_Include_IBF_1"));
    {
        load_ibf(ibf, "IBF_Test_Include_IBF_1");
        auto agent = ibf.membership_agent();

        std::vector<bool> expected_result(1, 0);
        auto & res = agent.bulk_contains(2);
        EXPECT_RANGE_EQ(expected_result,  res);
        expected_result[0] = 1;
        auto & res2 = agent.bulk_contains(24);
        EXPECT_RANGE_EQ(expected_result,  res2);
    }
}

TEST_F(ibf_test, given_expression_thresholds_exclude_file)
{
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    initialization_args(ibf_args);
    ibf_args.path_out = "IBF_Test_Exclude_";
    ibf_args.expression_thresholds = {1, 2};
    minimiser_args.exclude_file = data("mini_gen.fasta");
    std::vector<std::filesystem::path> sequence_files = {data("mini_example.fasta")};
    std::vector<double> fpr = {0.05};

    std::vector<uint16_t> expected{1, 2};
    std::vector<uint8_t> cutoffs{0};

    std::vector<uint16_t> medians = ibf(sequence_files, ibf_args, minimiser_args, fpr, cutoffs);

    EXPECT_EQ(expected, medians);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    ASSERT_TRUE(std::filesystem::exists("IBF_Test_Exclude_IBF_1"));
    {
        load_ibf(ibf, "IBF_Test_Exclude_IBF_1");
        auto agent = ibf.membership_agent();

        std::vector<bool> expected_result(1, 0);
        auto & res = agent.bulk_contains(2);
        EXPECT_RANGE_EQ(expected_result,  res);
        expected_result[0] = 1;
        auto & res2 = agent.bulk_contains(24);
        EXPECT_RANGE_EQ(expected_result,  res2);
    }
}

TEST_F(ibf_test, no_given_expression_thresholds)
{
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    initialization_args(ibf_args);
    ibf_args.path_out = "IBF_Test_";
    ibf_args.number_expression_thresholds = 2;
    std::vector<std::filesystem::path> sequence_files = {data("mini_example.fasta")};
    std::vector<double> fpr = {0.05};

    std::vector<uint16_t> expected{};
    std::vector<uint8_t> cutoffs{0};

    std::vector<uint16_t> medians = ibf(sequence_files, ibf_args, minimiser_args, fpr, cutoffs);

    EXPECT_EQ(expected, medians);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;

    ASSERT_TRUE(std::filesystem::exists("IBF_Test_IBF_Level_0"));
    {
        load_ibf(ibf, "IBF_Test_IBF_Level_0");
        auto agent = ibf.membership_agent();

        std::vector<bool> expected_result(1, 0);
        auto & res = agent.bulk_contains(2);
        EXPECT_RANGE_EQ(expected_result,  res);
        expected_result[0] = 1;
        auto & res2 = agent.bulk_contains(24);
        EXPECT_RANGE_EQ(expected_result,  res2);
    }
}

TEST_F(ibf_test, expression_thresholds_by_genome)
{
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    initialization_args(ibf_args);
    ibf_args.path_out = "IBF_Test_";
    ibf_args.number_expression_thresholds = 1;
    std::vector<std::filesystem::path> sequence_files = {data("mini_example.fasta")};
    std::vector<double> fpr = {0.05};

    std::vector<uint16_t> expected{};
    std::vector<uint8_t> cutoffs{};

    std::vector<uint16_t> medians = ibf(sequence_files, ibf_args, minimiser_args, fpr, cutoffs,
                                        data("mini_gen.fasta"));

    EXPECT_EQ(expected, medians);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;

    ASSERT_TRUE(std::filesystem::exists("IBF_Test_IBF_Level_0"));
    {
        load_ibf(ibf, "IBF_Test_IBF_Level_0");
        auto agent = ibf.membership_agent();

        std::vector<bool> expected_result(1, 0);
        auto & res = agent.bulk_contains(2);
        EXPECT_RANGE_EQ(expected_result,  res);
        expected_result[0] = 1;
        auto & res2 = agent.bulk_contains(192);
        EXPECT_RANGE_EQ(expected_result,  res2);
    }
}

TEST_F(ibf_test, throws)
{
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    std::vector<uint8_t> cutoffs{};
    initialization_args(ibf_args);
    ibf_args.path_out = "IBF_Test_";
    std::vector<std::filesystem::path> sequence_files = {data("mini_example.fasta")};
    std::vector<double> fpr = {0.05};

    EXPECT_THROW(ibf(sequence_files, ibf_args, minimiser_args, fpr, cutoffs), std::invalid_argument);

    ibf_args.number_expression_thresholds = 0;
    fpr = {};
    EXPECT_THROW(ibf(sequence_files, ibf_args, minimiser_args, fpr, cutoffs), std::invalid_argument);

    fpr = {0.05};
    EXPECT_THROW(ibf(sequence_files, ibf_args, minimiser_args, fpr, cutoffs), std::invalid_argument);

    ibf_args.expression_thresholds = {10, 1000};
    EXPECT_THROW(ibf(sequence_files, ibf_args, minimiser_args, fpr, cutoffs), std::invalid_argument);
}

TEST_F(ibf_test, given_cutoffs)
{
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    initialization_args(ibf_args);
    ibf_args.path_out = "IBF_Test_Cut_";
    ibf_args.expression_thresholds = {1, 2};
    std::vector<uint8_t> cutoffs = {0};
    std::vector<std::filesystem::path> sequence_files = {data("mini_example.fasta")};
    std::vector<double> fpr = {0.05};

    std::vector<uint16_t> expected{1, 2};

    std::vector<uint16_t> medians = ibf(sequence_files, ibf_args, minimiser_args, fpr, cutoffs);

    EXPECT_EQ(expected, medians);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    ASSERT_TRUE(std::filesystem::exists("IBF_Test_Cut_IBF_1"));
    {
        load_ibf(ibf, "IBF_Test_Cut_IBF_1");
        auto agent = ibf.membership_agent();

        std::vector<bool> expected_result(1, 0);
        auto & res = agent.bulk_contains(2);
        EXPECT_RANGE_EQ(expected_result,  res);
        expected_result[0] = 1;
        auto & res2 = agent.bulk_contains(24);
        EXPECT_RANGE_EQ(expected_result,  res2);
    }

    estimate_ibf_arguments args{};
    ASSERT_TRUE(std::filesystem::exists("IBF_Test_Cut_IBF_Data"));
    {
        load_args(args, "IBF_Test_Cut_IBF_Data");
        EXPECT_EQ(4, args.k);
        EXPECT_EQ(4, args.w_size.get());
        EXPECT_EQ(0, args.s.get());
        EXPECT_EQ(15, args.shape.to_ulong());
        EXPECT_EQ(true, args.compressed);
        EXPECT_EQ(ibf_args.number_expression_thresholds, args.number_expression_thresholds);
        EXPECT_RANGE_EQ(ibf_args.expression_thresholds, args.expression_thresholds);
        EXPECT_EQ(ibf_args.samplewise, args.samplewise);
    }
}

TEST_F(ibf_test, different_file_sizes)
{
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    initialization_args(ibf_args);
    ibf_args.path_out = "IBF_Test_Diff_";
    ibf_args.number_expression_thresholds = 4;
    std::vector<std::filesystem::path> sequence_files = {data("mini_example.fasta"), data("exp_01.fasta")};
    std::vector<double> fpr = {0.05};

    std::vector<uint16_t> expected{};
    std::vector<uint8_t> cutoffs{};

    std::vector<uint16_t> medians = ibf(sequence_files, ibf_args, minimiser_args, fpr, cutoffs);

    EXPECT_EQ(expected, medians);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;

    ASSERT_TRUE(std::filesystem::exists("IBF_Test_Diff_IBF_Level_0"));
    {
        load_ibf(ibf, "IBF_Test_Diff_IBF_Level_0");
        auto agent = ibf.membership_agent();

        std::vector<bool> expected_result(2, 0);
        auto & res = agent.bulk_contains(2);
        EXPECT_RANGE_EQ(expected_result,  res);
        expected_result[0] = 1;
        expected_result[1] = 1;
        auto & res2 = agent.bulk_contains(24);
        EXPECT_RANGE_EQ(expected_result,  res2);
    }
}
