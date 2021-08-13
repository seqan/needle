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

std::vector<robin_hood::unordered_node_map<uint64_t, uint16_t>> expected_hash_tables{   // minimisers:
                                                                             {{0,2},   // AAAA
                                                                              {1,4},   // AAAC
                                                                              {6,4},   // AACG
                                                                              {24,1},  // ACGA
                                                                              {27,5},  // ACGT
                                                                              {97,3},  // CGAC
                                                                              {108,2}, // CGTA
                                                                              {109,3}, // CGTC
                                                                              {112,3}, // CTAA
                                                                              {177,1}, // GTAC
                                                                              {192,3}, // TAAA
                                                                              {216,1}, // TCGA
                                                                                     },
                                                                             {{27,1},  // ACGT
                                                                              {42,1},  // AGGG
                                                                              {74,1},  // CAGG
                                                                              {82,1},  // CCAG
                                                                              {84,1},  // CCCA
                                                                              {85,19}, // CCCC
                                                                              {86,1},  // CCCG
                                                                              {109,1}, // CGTC
                                                                              {149,2}, // GCCC
                                                                              {161,1}, // GGAC
                                                                              {165,1}, // GGCC
                                                                              {168,1}, // GGGA
                                                                                     },};

void initialization_args(estimate_ibf_arguments & args)
{
    args.k = 4;
    args.shape = seqan3::ungapped{args.k};
    args.w_size = seqan3::window_size{4};
    args.s = seqan3::seed{0};
    args.path_out = tmp_dir/"Test_";
    args.fpr = {0.05};
    args.compressed = true;
}

TEST(minimiser, small_example)
{
    estimate_ibf_arguments args{};
    minimiser_arguments minimiser_args{};
    initialization_args(args);
    minimiser_args.cutoffs = {0, 0};
    args.expression_levels = {0};
    std::vector<std::filesystem::path> sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta",
                                                         std::string(DATA_INPUT_DIR) + "mini_example2.fasta"};
    minimiser(sequence_files, args, minimiser_args);
    uint32_t normalized_exp_value{};
    robin_hood::unordered_node_map<uint64_t, uint16_t> result_hash_table{};
    std::vector<std::filesystem::path> minimiser_files{};
    uint64_t num_of_minimisers{};
    std::vector<uint64_t> expected_nums{12, 12};

    for (int i = 0; i < sequence_files.size(); ++i)
    {
        // Test Header file
        read_binary_start(args, tmp_dir/("Test_" + std::string{sequence_files[i].stem()} + ".minimiser"), num_of_minimisers);

        EXPECT_EQ(4, args.k);
        EXPECT_EQ(4, args.w_size.get());
        EXPECT_EQ(0, args.s.get());
        EXPECT_EQ(15, args.shape.to_ulong());
        EXPECT_EQ(expected_nums[i], num_of_minimisers);

        // Test binary file
        read_binary(tmp_dir/("Test_" + std::string{sequence_files[i].stem()} + ".minimiser"), result_hash_table);
        minimiser_files.push_back(tmp_dir/("Test_" + std::string{sequence_files[i].stem()} + ".minimiser"));
        for (auto & hash : expected_hash_tables[i])
        {
            EXPECT_EQ(expected_hash_tables[i][hash.first], result_hash_table[hash.first]);
        }

        result_hash_table.clear();
    }

    EXPECT_EQ(args.expression_levels, ibf(minimiser_files, args));
    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    load_ibf(ibf, tmp_dir/"Test_IBF_0");
    auto agent = ibf.membership_agent();

    std::vector<bool> expected_result(2, 0);
    auto & res = agent.bulk_contains(2);
    EXPECT_RANGE_EQ(expected_result,  res);
    expected_result[0] = 1;
    auto & res2 = agent.bulk_contains(0);
    EXPECT_RANGE_EQ(expected_result,  res2);
    expected_result[1] = 1;
    auto & res3 = agent.bulk_contains(27);
    EXPECT_RANGE_EQ(expected_result,  res3);

    std::filesystem::remove(tmp_dir/"Test_IBF_0");
    std::filesystem::remove(tmp_dir/("Test_mini_example.minimiser"));
    std::filesystem::remove(tmp_dir/("Test_mini_example2.minimiser"));
}

TEST(minimiser, small_example_different_shape)
{
    estimate_ibf_arguments args{};
    minimiser_arguments minimiser_args{};
    initialization_args(args);
    minimiser_args.cutoffs = {0, 0};
    args.shape = seqan3::bin_literal{14};
    std::vector<std::filesystem::path> sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta",
                                                         std::string(DATA_INPUT_DIR) + "mini_example2.fasta"};
    minimiser(sequence_files, args, minimiser_args);

    uint64_t num_of_minimisers{};
    std::vector<uint64_t> expected_nums{8, 10};

    for (int i = 0; i < sequence_files.size(); ++i)
    {
        // Test Header file
        read_binary_start(args, tmp_dir/("Test_" + std::string{sequence_files[i].stem()} + ".minimiser"), num_of_minimisers);

        EXPECT_EQ(4, args.k);
        EXPECT_EQ(4, args.w_size.get());
        EXPECT_EQ(0, args.s.get());
        EXPECT_EQ(14, args.shape.to_ulong());
        EXPECT_EQ(expected_nums[i], num_of_minimisers);
    }

    std::filesystem::remove(tmp_dir/("Test_mini_example.minimiser"));
    std::filesystem::remove(tmp_dir/("Test_mini_example2.minimiser"));
}

TEST(minimiser, small_example_samplewise)
{
    estimate_ibf_arguments args{};
    minimiser_arguments minimiser_args{};
    initialization_args(args);

    minimiser_args.cutoffs = {0, 0};
    args.number_expression_levels = 1;
    std::vector<std::filesystem::path> sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta",
                                                         std::string(DATA_INPUT_DIR) + "mini_example2.fasta"};

    minimiser(sequence_files, args, minimiser_args);
    uint32_t normalized_exp_value{};
    std::vector<std::vector<uint32_t>> expected_counts{{7}, {12}};
    std::vector<uint16_t> expected_levels{};
    robin_hood::unordered_node_map<uint64_t, uint16_t> result_hash_table{};
    std::vector<std::filesystem::path> minimiser_files{};
    uint64_t num_of_minimisers{};
    std::vector<uint64_t> expected_nums{12, 12};

    for (int i = 0; i < sequence_files.size(); ++i)
    {
        // Test Header file
        args.expression_levels = {};
        read_binary_start(args, tmp_dir/("Test_" + std::string{sequence_files[i].stem()} + ".minimiser"), num_of_minimisers);
        EXPECT_EQ(4, args.k);
        EXPECT_EQ(4, args.w_size.get());
        EXPECT_EQ(0, args.s.get());
        EXPECT_EQ(15, args.shape.to_ulong());
        EXPECT_EQ(expected_nums[i], num_of_minimisers);

        // Test binary file
        read_binary(tmp_dir/("Test_" + std::string{sequence_files[i].stem()} + ".minimiser"), result_hash_table);
        minimiser_files.push_back(tmp_dir/("Test_" + std::string{sequence_files[i].stem()} + ".minimiser"));
        for (auto & hash : expected_hash_tables[i])
        {
            EXPECT_EQ(expected_hash_tables[i][hash.first], result_hash_table[hash.first]);
        }

        result_hash_table.clear();
    }
    args.expression_levels = {};
    EXPECT_EQ(expected_levels, ibf(minimiser_files, args));

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    load_ibf(ibf, tmp_dir/"Test_IBF_Level_0");
    auto agent = ibf.membership_agent();

    std::vector<bool> expected_result(2, 0);
    auto & res = agent.bulk_contains(2);
    EXPECT_RANGE_EQ(expected_result,  res);
    auto & res2 = agent.bulk_contains(0);
    EXPECT_RANGE_EQ(expected_result,  res2);
    expected_result[0] = 1;
    expected_result[1] = 1;
    auto & res3 = agent.bulk_contains(27);
    EXPECT_RANGE_EQ(expected_result,  res3);
    std::filesystem::remove(tmp_dir/"Test_IBF_Level_0");
    std::filesystem::remove(tmp_dir/"Test_IBF_Levels.levels");
    std::filesystem::remove(tmp_dir/("Test_mini_example.minimiser"));
    std::filesystem::remove(tmp_dir/("Test_mini_example2.minimiser"));
}

TEST(minimiser, cutoff_by_filesize)
{
    estimate_ibf_arguments args{};
    minimiser_arguments minimiser_args{};
    initialization_args(args);
    args.expression_levels = {0};
    std::vector<std::filesystem::path> sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta",
                                                         std::string(DATA_INPUT_DIR) + "mini_example2.fasta"};

    minimiser(sequence_files, args, minimiser_args);
    uint32_t normalized_exp_value{};
    std::vector<std::filesystem::path> minimiser_files{};
    uint64_t num_of_minimisers{};
    std::vector<uint64_t> expected_nums{12, 12};

    for (int i = 0; i < sequence_files.size(); ++i)
    {
        // Test Header file
        read_binary_start(args, tmp_dir/("Test_" + std::string{sequence_files[i].stem()} + ".minimiser"), num_of_minimisers);

        EXPECT_EQ(4, args.k);
        EXPECT_EQ(4, args.w_size.get());
        EXPECT_EQ(0, args.s.get());
        EXPECT_EQ(15, args.shape.to_ulong());
        EXPECT_EQ(expected_nums[i], num_of_minimisers);
        minimiser_files.push_back(tmp_dir/("Test_" + std::string{sequence_files[i].stem()} + ".minimiser"));
    }

    EXPECT_EQ(args.expression_levels, ibf(minimiser_files, args));

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    load_ibf(ibf, tmp_dir/"Test_IBF_0");
    auto agent = ibf.membership_agent();

    std::vector<bool> expected_result(2, 0);
    auto & res = agent.bulk_contains(2);
    EXPECT_RANGE_EQ(expected_result,  res);
    expected_result[0] = 1;
    auto & res2 = agent.bulk_contains(0);
    EXPECT_RANGE_EQ(expected_result, res2);
    expected_result[0] = 0;
    expected_result[1] = 1;
    auto & res3 = agent.bulk_contains(85);
    EXPECT_RANGE_EQ(expected_result, res3);

    std::filesystem::remove(tmp_dir/"Test_IBF_0");
    std::filesystem::remove(tmp_dir/("Test_mini_example.minimiser"));
    std::filesystem::remove(tmp_dir/("Test_mini_example2.minimiser"));
}

TEST(minimiser, small_example_two_threads)
{
    estimate_ibf_arguments args{};
    minimiser_arguments minimiser_args{};
    initialization_args(args);
    args.threads = 2;
    minimiser_args.cutoffs = {0, 0};
    args.expression_levels = {0};
    std::vector<std::filesystem::path> sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta",
                                                         std::string(DATA_INPUT_DIR) + "mini_example2.fasta"};
    minimiser(sequence_files, args, minimiser_args);
    args.threads = 1;
    uint32_t normalized_exp_value{};
    robin_hood::unordered_node_map<uint64_t, uint16_t> result_hash_table{};
    std::vector<std::filesystem::path> minimiser_files{};
    uint64_t num_of_minimisers{};
    std::vector<uint64_t> expected_nums{12, 12};

    for (int i = 0; i < sequence_files.size(); ++i)
    {
        // Test Header file
        read_binary_start(args, tmp_dir/("Test_" + std::string{sequence_files[i].stem()} + ".minimiser"), num_of_minimisers);

        EXPECT_EQ(4, args.k);
        EXPECT_EQ(4, args.w_size.get());
        EXPECT_EQ(0, args.s.get());
        EXPECT_EQ(15, args.shape.to_ulong());
        EXPECT_EQ(expected_nums[i], num_of_minimisers);

        // Test binary file
        read_binary(tmp_dir/("Test_" + std::string{sequence_files[i].stem()} + ".minimiser"), result_hash_table);
        minimiser_files.push_back(tmp_dir/("Test_" + std::string{sequence_files[i].stem()} + ".minimiser"));
        for (auto & hash : expected_hash_tables[i])
        {
            EXPECT_EQ(expected_hash_tables[i][hash.first], result_hash_table[hash.first]);
        }

        result_hash_table.clear();
    }

    EXPECT_EQ(args.expression_levels, ibf(minimiser_files, args));

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    load_ibf(ibf, tmp_dir/"Test_IBF_0");
    auto agent = ibf.membership_agent();

    std::vector<bool> expected_result(2, 0);
    auto & res = agent.bulk_contains(2);
    EXPECT_RANGE_EQ(expected_result,  res);
    expected_result[0] = 1;
    auto & res2 = agent.bulk_contains(0);
    EXPECT_RANGE_EQ(expected_result, res2);
    expected_result[1] = 1;
    auto & res3 = agent.bulk_contains(27);
    EXPECT_RANGE_EQ(expected_result, res3);

    std::filesystem::remove(tmp_dir/"Test_IBF_0");
    std::filesystem::remove(tmp_dir/("Test_mini_example.minimiser"));
    std::filesystem::remove(tmp_dir/("Test_mini_example2.minimiser"));
}
