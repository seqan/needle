// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <iostream>

#include <seqan3/test/expect_range_eq.hpp>

#include "../app_test.hpp"
#include "ibf.hpp"
#include "minimiser.hpp"
#include "misc/stream.hpp"
#include "shared.hpp"

// To prevent issues when running multiple CLI tests in parallel, give each CLI test unique names:
struct minimiser_test : public app_test
{
    void initialization_args(estimate_ibf_arguments & args)
    {
        args.k = 4;
        args.shape = seqan3::ungapped{args.k};
        args.w_size = seqan3::window_size{4};
        args.s = seqan3::seed{0};
        args.path_out = "Minimiser_Test_";
    }

    std::vector<robin_hood::unordered_node_map<uint64_t, uint16_t>> expected_hash_tables{
        // minimisers:
        {
            {0, 2},   // AAAA
            {1, 4},   // AAAC
            {6, 4},   // AACG
            {24, 1},  // ACGA
            {27, 5},  // ACGT
            {97, 3},  // CGAC
            {108, 2}, // CGTA
            {109, 3}, // CGTC
            {112, 3}, // CTAA
            {177, 1}, // GTAC
            {192, 3}, // TAAA
            {216, 1}, // TCGA
        },
        {
            {27, 1},  // ACGT
            {42, 1},  // AGGG
            {74, 1},  // CAGG
            {82, 1},  // CCAG
            {84, 1},  // CCCA
            {85, 19}, // CCCC
            {86, 1},  // CCCG
            {109, 1}, // CGTC
            {149, 2}, // GCCC
            {161, 1}, // GGAC
            {165, 1}, // GGCC
            {168, 1}, // GGGA
        },
    };
};

TEST_F(minimiser_test, small_example)
{
    estimate_ibf_arguments args{};
    minimiser_file_input_arguments minimiser_args{};
    initialization_args(args);
    std::vector<uint8_t> cutoffs = {0, 0};
    args.expression_thresholds = {0};
    std::vector<double> fpr = {0.05};
    std::vector<std::filesystem::path> sequence_files = {data("mini_example.fasta"), data("mini_example2.fasta")};
    minimiser(sequence_files, args, minimiser_args, cutoffs);
    robin_hood::unordered_node_map<uint64_t, uint16_t> result_hash_table{};
    std::vector<std::filesystem::path> minimiser_files{};
    uint64_t num_of_minimisers{};
    std::vector<uint64_t> expected_nums{12, 12};

    for (size_t i = 0; i < sequence_files.size(); ++i)
    {
        uint8_t cutoff{};
        // Test Header file
        read_binary_start(args,
                          ("Minimiser_Test_" + std::string{sequence_files[i].stem()} + ".minimiser"),
                          num_of_minimisers,
                          cutoff);

        EXPECT_EQ(4, args.k);
        EXPECT_EQ(4, args.w_size.get());
        EXPECT_EQ(0, args.s.get());
        EXPECT_EQ(15, args.shape.to_ulong());
        EXPECT_EQ(expected_nums[i], num_of_minimisers);
        EXPECT_EQ(0, cutoff);

        // Test binary file
        read_binary(("Minimiser_Test_" + std::string{sequence_files[i].stem()} + ".minimiser"), result_hash_table);
        minimiser_files.push_back(("Minimiser_Test_" + std::string{sequence_files[i].stem()} + ".minimiser"));
        for (auto & hash : expected_hash_tables[i])
        {
            EXPECT_EQ(expected_hash_tables[i][hash.first], result_hash_table[hash.first]);
        }

        result_hash_table.clear();
    }

    EXPECT_EQ(args.expression_thresholds, ibf(minimiser_files, args, fpr));
    seqan::hibf::interleaved_bloom_filter ibf;
    load_ibf(ibf, "Minimiser_Test_IBF_0");
    auto agent = ibf.containment_agent();

    std::vector<bool> expected_result(2, 0);
    auto & res = agent.bulk_contains(2);
    EXPECT_RANGE_EQ(expected_result, res);
    expected_result[0] = 1;
    auto & res2 = agent.bulk_contains(0);
    EXPECT_RANGE_EQ(expected_result, res2);
    expected_result[1] = 1;
    auto & res3 = agent.bulk_contains(27);
    EXPECT_RANGE_EQ(expected_result, res3);
}

TEST_F(minimiser_test, small_example_different_shape)
{
    estimate_ibf_arguments args{};
    minimiser_file_input_arguments minimiser_args{};
    initialization_args(args);
    std::vector<uint8_t> cutoffs = {0, 0};
    args.shape = seqan3::bin_literal{0b1101};
    EXPECT_EQ(13, args.shape.to_ulong());
    std::vector<std::filesystem::path> sequence_files = {data("mini_example.fasta"), data("mini_example2.fasta")};
    minimiser(sequence_files, args, minimiser_args, cutoffs);

    uint64_t num_of_minimisers{};
    std::vector<uint64_t> expected_nums{12, 11};

    for (size_t i = 0; i < sequence_files.size(); ++i)
    {
        uint8_t cutoff{};
        // Test Header file
        read_binary_start(args,
                          ("Minimiser_Test_" + std::string{sequence_files[i].stem()} + ".minimiser"),
                          num_of_minimisers,
                          cutoff);

        EXPECT_EQ(4, args.k);
        EXPECT_EQ(4, args.w_size.get());
        EXPECT_EQ(0, args.s.get());
        EXPECT_EQ(13, args.shape.to_ulong());
        EXPECT_EQ(expected_nums[i], num_of_minimisers);
        EXPECT_EQ(0, cutoff);
    }
}

TEST_F(minimiser_test, small_example_samplewise)
{
    estimate_ibf_arguments args{};
    minimiser_file_input_arguments minimiser_args{};
    initialization_args(args);

    std::vector<uint8_t> cutoffs = {0, 0};
    args.number_expression_thresholds = 1;
    std::vector<double> fpr = {0.05};
    std::vector<std::filesystem::path> sequence_files = {data("mini_example.fasta"), data("mini_example2.fasta")};

    minimiser(sequence_files, args, minimiser_args, cutoffs);
    std::vector<std::vector<uint32_t>> expected_counts{{7}, {12}};
    std::vector<uint16_t> expected_levels{};
    robin_hood::unordered_node_map<uint64_t, uint16_t> result_hash_table{};
    std::vector<std::filesystem::path> minimiser_files{};
    uint64_t num_of_minimisers{};
    std::vector<uint64_t> expected_nums{12, 12};

    for (size_t i = 0; i < sequence_files.size(); ++i)
    {
        uint8_t cutoff{};
        // Test Header file
        args.expression_thresholds = {};
        read_binary_start(args,
                          ("Minimiser_Test_" + std::string{sequence_files[i].stem()} + ".minimiser"),
                          num_of_minimisers,
                          cutoff);
        EXPECT_EQ(4, args.k);
        EXPECT_EQ(4, args.w_size.get());
        EXPECT_EQ(0, args.s.get());
        EXPECT_EQ(15, args.shape.to_ulong());
        EXPECT_EQ(expected_nums[i], num_of_minimisers);
        EXPECT_EQ(0, cutoff);

        // Test binary file
        read_binary(("Minimiser_Test_" + std::string{sequence_files[i].stem()} + ".minimiser"), result_hash_table);
        minimiser_files.push_back(("Minimiser_Test_" + std::string{sequence_files[i].stem()} + ".minimiser"));
        for (auto & hash : expected_hash_tables[i])
        {
            EXPECT_EQ(expected_hash_tables[i][hash.first], result_hash_table[hash.first]);
        }

        result_hash_table.clear();
    }
    args.expression_thresholds = {};
    EXPECT_EQ(expected_levels, ibf(minimiser_files, args, fpr));

    seqan::hibf::interleaved_bloom_filter ibf;
    load_ibf(ibf, "Minimiser_Test_IBF_Level_0");
    auto agent = ibf.containment_agent();

    std::vector<bool> expected_result(2, 0);
    auto & res = agent.bulk_contains(2);
    EXPECT_RANGE_EQ(expected_result, res);
    auto & res2 = agent.bulk_contains(0);
    expected_result[0] = 1;
    EXPECT_RANGE_EQ(expected_result, res2);
    expected_result[1] = 1;
    auto & res3 = agent.bulk_contains(27);
    EXPECT_RANGE_EQ(expected_result, res3);
}

TEST_F(minimiser_test, cutoff_by_filesize)
{
    estimate_ibf_arguments args{};
    minimiser_file_input_arguments minimiser_args{};
    initialization_args(args);
    args.expression_thresholds = {0};
    std::vector<double> fpr = {0.05};
    std::vector<uint8_t> cutoffs{};
    std::vector<std::filesystem::path> sequence_files = {data("mini_example.fasta"), data("mini_example2.fasta")};

    minimiser(sequence_files, args, minimiser_args, cutoffs);

    std::vector<std::filesystem::path> minimiser_files{};
    uint64_t num_of_minimisers{};
    std::vector<uint64_t> expected_nums{12, 12};

    for (size_t i = 0; i < sequence_files.size(); ++i)
    {
        uint8_t cutoff{};
        // Test Header file
        read_binary_start(args,
                          ("Minimiser_Test_" + std::string{sequence_files[i].stem()} + ".minimiser"),
                          num_of_minimisers,
                          cutoff);

        EXPECT_EQ(4, args.k);
        EXPECT_EQ(4, args.w_size.get());
        EXPECT_EQ(0, args.s.get());
        EXPECT_EQ(15, args.shape.to_ulong());
        EXPECT_EQ(0, cutoff);
        EXPECT_EQ(expected_nums[i], num_of_minimisers);
        minimiser_files.push_back(("Minimiser_Test_" + std::string{sequence_files[i].stem()} + ".minimiser"));
    }

    EXPECT_EQ(args.expression_thresholds, ibf(minimiser_files, args, fpr));

    seqan::hibf::interleaved_bloom_filter ibf;
    load_ibf(ibf, "Minimiser_Test_IBF_0");
    auto agent = ibf.containment_agent();

    std::vector<bool> expected_result(2, 0);
    auto & res = agent.bulk_contains(2);
    EXPECT_RANGE_EQ(expected_result, res);
    expected_result[0] = 1;
    auto & res2 = agent.bulk_contains(0);
    EXPECT_RANGE_EQ(expected_result, res2);
    expected_result[0] = 0;
    expected_result[1] = 1;
    auto & res3 = agent.bulk_contains(85);
    EXPECT_RANGE_EQ(expected_result, res3);
}

TEST_F(minimiser_test, small_example_two_threads)
{
    estimate_ibf_arguments args{};
    minimiser_file_input_arguments minimiser_args{};
    initialization_args(args);
    args.threads = 2;
    std::vector<uint8_t> cutoffs = {0, 0};
    args.expression_thresholds = {0};
    std::vector<double> fpr = {0.05};
    std::vector<std::filesystem::path> sequence_files = {data("mini_example.fasta"), data("mini_example2.fasta")};
    minimiser(sequence_files, args, minimiser_args, cutoffs);
    args.threads = 1;
    robin_hood::unordered_node_map<uint64_t, uint16_t> result_hash_table{};
    std::vector<std::filesystem::path> minimiser_files{};
    uint64_t num_of_minimisers{};
    std::vector<uint64_t> expected_nums{12, 12};

    for (size_t i = 0; i < sequence_files.size(); ++i)
    {
        uint8_t cutoff{};
        // Test Header file
        read_binary_start(args,
                          ("Minimiser_Test_" + std::string{sequence_files[i].stem()} + ".minimiser"),
                          num_of_minimisers,
                          cutoff);

        EXPECT_EQ(4, args.k);
        EXPECT_EQ(4, args.w_size.get());
        EXPECT_EQ(0, args.s.get());
        EXPECT_EQ(15, args.shape.to_ulong());
        EXPECT_EQ(0, cutoff);
        EXPECT_EQ(expected_nums[i], num_of_minimisers);

        // Test binary file
        read_binary(("Minimiser_Test_" + std::string{sequence_files[i].stem()} + ".minimiser"), result_hash_table);
        minimiser_files.push_back(("Minimiser_Test_" + std::string{sequence_files[i].stem()} + ".minimiser"));
        for (auto & hash : expected_hash_tables[i])
        {
            EXPECT_EQ(expected_hash_tables[i][hash.first], result_hash_table[hash.first]);
        }

        result_hash_table.clear();
    }

    EXPECT_EQ(args.expression_thresholds, ibf(minimiser_files, args, fpr));

    seqan::hibf::interleaved_bloom_filter ibf;
    load_ibf(ibf, "Minimiser_Test_IBF_0");
    auto agent = ibf.containment_agent();

    std::vector<bool> expected_result(2, 0);
    auto & res = agent.bulk_contains(2);
    EXPECT_RANGE_EQ(expected_result, res);
    expected_result[0] = 1;
    auto & res2 = agent.bulk_contains(0);
    EXPECT_RANGE_EQ(expected_result, res2);
    expected_result[1] = 1;
    auto & res3 = agent.bulk_contains(27);
    EXPECT_RANGE_EQ(expected_result, res3);
}

TEST_F(minimiser_test, small_example_include)
{
    estimate_ibf_arguments args{};
    minimiser_file_input_arguments minimiser_args{};
    initialization_args(args);
    args.path_out = "Minimiser_Test_In_";
    std::vector<uint8_t> cutoffs = {0, 0};
    minimiser_args.include_file = data("mini_gen.fasta");
    std::vector<std::filesystem::path> sequence_files = {data("mini_example.fasta"), data("mini_example2.fasta")};
    minimiser(sequence_files, args, minimiser_args, cutoffs);
    robin_hood::unordered_node_map<uint64_t, uint16_t> result_hash_table{};
    std::vector<std::filesystem::path> minimiser_files{};
    uint64_t num_of_minimisers{};
    std::vector<uint64_t> expected_nums{1, 0};

    for (size_t i = 0; i < sequence_files.size(); ++i)
    {
        uint8_t cutoff{};
        // Test Header file
        read_binary_start(args,
                          ("Minimiser_Test_In_" + std::string{sequence_files[i].stem()} + ".minimiser"),
                          num_of_minimisers,
                          cutoff);

        EXPECT_EQ(4, args.k);
        EXPECT_EQ(4, args.w_size.get());
        EXPECT_EQ(0, args.s.get());
        EXPECT_EQ(15, args.shape.to_ulong());
        EXPECT_EQ(0, cutoff);
        EXPECT_EQ(expected_nums[i], num_of_minimisers);

        // Test binary file
        read_binary(("Minimiser_Test_In_" + std::string{sequence_files[i].stem()} + ".minimiser"), result_hash_table);
        minimiser_files.push_back(("Minimiser_Test_In_" + std::string{sequence_files[i].stem()} + ".minimiser"));
        if (i == 0)
        {
            for (auto & hash : result_hash_table)
            {
                EXPECT_EQ(192, hash.first); // 192 minimiser TAAA, only minimiser in mini_gen
                EXPECT_EQ(3, result_hash_table[hash.first]);
            }
        }
        else
        {
            EXPECT_EQ(0, result_hash_table.size());
        }

        result_hash_table.clear();
    }
}

TEST_F(minimiser_test, small_example_exclude)
{
    estimate_ibf_arguments args{};
    minimiser_file_input_arguments minimiser_args{};
    initialization_args(args);
    args.path_out = "Minimiser_Test_Ex_";
    std::vector<uint8_t> cutoffs = {0, 0};
    minimiser_args.exclude_file = data("mini_gen2.fasta");
    args.expression_thresholds = {0};
    std::vector<double> fpr = {0.05};
    std::vector<std::filesystem::path> sequence_files = {data("mini_example.fasta"), data("mini_example2.fasta")};
    minimiser(sequence_files, args, minimiser_args, cutoffs);
    robin_hood::unordered_node_map<uint64_t, uint16_t> result_hash_table{};
    std::vector<std::filesystem::path> minimiser_files{};
    uint64_t num_of_minimisers{};
    std::vector<uint64_t> expected_nums{11, 12};

    for (size_t i = 0; i < sequence_files.size(); ++i)
    {
        uint8_t cutoff{};
        // Test Header file
        read_binary_start(args,
                          ("Minimiser_Test_Ex_" + std::string{sequence_files[i].stem()} + ".minimiser"),
                          num_of_minimisers,
                          cutoff);

        EXPECT_EQ(4, args.k);
        EXPECT_EQ(4, args.w_size.get());
        EXPECT_EQ(0, args.s.get());
        EXPECT_EQ(15, args.shape.to_ulong());
        EXPECT_EQ(expected_nums[i], num_of_minimisers);
        EXPECT_EQ(0, cutoff);

        // Test binary file
        read_binary(("Minimiser_Test_Ex_" + std::string{sequence_files[i].stem()} + ".minimiser"), result_hash_table);
        minimiser_files.push_back(("Minimiser_Test_Ex_" + std::string{sequence_files[i].stem()} + ".minimiser"));
        EXPECT_FALSE(result_hash_table.contains(216)); //TCGA, only minimiser in mini_gen2
        for (auto & hash : expected_hash_tables[i])
        {
            if (hash.first == 216) //TCGA, only minimiser in mini_gen2
                continue;

            EXPECT_EQ(expected_hash_tables[i][hash.first], result_hash_table[hash.first]);
        }

        result_hash_table.clear();
    }

    EXPECT_EQ(args.expression_thresholds, ibf(minimiser_files, args, fpr));
    seqan::hibf::interleaved_bloom_filter ibf;
    load_ibf(ibf, "Minimiser_Test_Ex_IBF_0");
    auto agent = ibf.containment_agent();

    std::vector<bool> expected_result(2, 0);
    auto & res = agent.bulk_contains(2);
    EXPECT_RANGE_EQ(expected_result, res);
    expected_result[0] = 1;
    auto & res2 = agent.bulk_contains(0);
    EXPECT_RANGE_EQ(expected_result, res2);
    expected_result[1] = 1;
    auto & res3 = agent.bulk_contains(27);
    EXPECT_RANGE_EQ(expected_result, res3);
}

TEST_F(minimiser_test, small_example_shape)
{
    std::vector<robin_hood::unordered_node_map<uint64_t, uint16_t>> expected_hash_tables_shape{
        // minimisers:
        {
            {0, 3},  // AA
            {1, 4},  // AC
            {2, 4},  // AG
            {3, 5},  // AT
            {4, 5},  // CA
            {5, 6},  // CC
            {9, 1},  // GC
            {12, 4}, // TA
        },
        {
            {2, 1},  // AT
            {3, 1},  // AG
            {4, 1},  // CA
            {5, 20}, // CC
            {6, 3},  // CG
            {8, 1},  // GA
            {9, 4},  // GC
        },
    };

    estimate_ibf_arguments args{};
    minimiser_file_input_arguments minimiser_args{};
    initialization_args(args);
    args.shape = seqan3::bin_literal{9};
    EXPECT_EQ(9, args.shape.to_ulong());
    args.path_out = "Minimiser_Test_Shape_";
    std::vector<uint8_t> cutoffs = {0, 0};
    args.expression_thresholds = {0};
    std::vector<double> fpr = {0.05};
    std::vector<std::filesystem::path> sequence_files = {data("mini_example.fasta"), data("mini_example2.fasta")};
    minimiser(sequence_files, args, minimiser_args, cutoffs);
    robin_hood::unordered_node_map<uint64_t, uint16_t> result_hash_table{};
    std::vector<std::filesystem::path> minimiser_files{};
    uint64_t num_of_minimisers{};
    std::vector<uint64_t> expected_nums{8, 7};

    for (size_t i = 0; i < sequence_files.size(); ++i)
    {
        uint8_t cutoff{};
        // Test Header file
        read_binary_start(args,
                          ("Minimiser_Test_Shape_" + std::string{sequence_files[i].stem()} + ".minimiser"),
                          num_of_minimisers,
                          cutoff);

        EXPECT_EQ(4, args.k);
        EXPECT_EQ(4, args.w_size.get());
        EXPECT_EQ(0, args.s.get());
        EXPECT_EQ(9, args.shape.to_ulong());
        EXPECT_EQ(expected_nums[i], num_of_minimisers);
        EXPECT_EQ(0, cutoff);

        // Test binary file
        read_binary(("Minimiser_Test_Shape_" + std::string{sequence_files[i].stem()} + ".minimiser"),
                    result_hash_table);
        minimiser_files.push_back(("Minimiser_Test_Shape_" + std::string{sequence_files[i].stem()} + ".minimiser"));
        for (auto & hash : expected_hash_tables_shape[i])
            EXPECT_EQ(expected_hash_tables_shape[i][hash.first], result_hash_table[hash.first]);

        result_hash_table.clear();
    }

    EXPECT_EQ(args.expression_thresholds, ibf(minimiser_files, args, fpr));
    seqan::hibf::interleaved_bloom_filter ibf;
    load_ibf(ibf, "Minimiser_Test_Shape_IBF_0");
    auto agent = ibf.containment_agent();

    std::vector<bool> expected_result(2, 0);
    auto & res = agent.bulk_contains(7);
    EXPECT_RANGE_EQ(expected_result, res);
    expected_result[0] = 1;
    auto & res2 = agent.bulk_contains(12);
    EXPECT_RANGE_EQ(expected_result, res2);
    expected_result[1] = 1;
    auto & res3 = agent.bulk_contains(2);
    EXPECT_RANGE_EQ(expected_result, res3);
}
