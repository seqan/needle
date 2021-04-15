#include <gtest/gtest.h>
#include <iostream>

#include "ibf.h"
#include "minimiser.h"
#include "estimate.h"

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

void initialization_args(arguments & args)
{
    args.compressed = true;
    args.k = 4;
    args.shape = seqan3::ungapped{args.k};
    args.w_size = seqan3::window_size{4};
    args.s = seqan3::seed{0};
}

void initialization_ibf_args(ibf_arguments & args)
{
    args.bin_size = {1000};
    args.path_out = tmp_dir/"Test_";
}

TEST(count, small_example)
{
    arguments args{};
    initialization_args(args);

    count(args, {std::string(DATA_INPUT_DIR) + "mini_example.fasta"}, std::string(DATA_INPUT_DIR) + "mini_gen.fasta",
          std::string(DATA_INPUT_DIR), false);

    std::ifstream output_file(tmp_dir/"mini_example.count.out");
    std::string line;
    std::string expected{"gen1\t3"};
    if (output_file.is_open())
    {
        while ( std::getline (output_file,line) )
        {
            EXPECT_EQ(expected,line);
        }
        output_file.close();
    }
    std::filesystem::remove(tmp_dir/"mini_example.count.out");
}

TEST(count, small_example_paired)
{
    arguments args{};
    initialization_args(args);

    count(args, {std::string(DATA_INPUT_DIR) + "mini_example.fasta", std::string(DATA_INPUT_DIR) + "mini_example.fasta"},
          std::string(DATA_INPUT_DIR) + "mini_gen.fasta",
          std::string(DATA_INPUT_DIR), true);

    std::ifstream output_file(tmp_dir/"mini_example.count.out");
    std::string line;
    std::string expected{"gen1\t6"};
    if (output_file.is_open())
    {
        while ( std::getline (output_file,line) )
        {
            EXPECT_EQ(expected,line);
        }
        output_file.close();
    }
    std::filesystem::remove(tmp_dir/"mini_example.count.out");
}

TEST(ibf, given_expression_levels)
{
    arguments args{};
    ibf_arguments ibf_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    ibf_args.expression_levels = {1, 2};
    ibf_args.sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};

    std::vector<uint32_t> expected{1, 2};

    std::vector<uint32_t> medians = ibf(args, ibf_args);

    EXPECT_EQ(expected, medians);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;

    load_ibf(ibf, tmp_dir/"Test_IBF_1");
    auto agent = ibf.membership_agent();

    sdsl::bit_vector expected_result(1, 0);
    EXPECT_EQ(expected_result,  agent.bulk_contains(2));
    expected_result[0] = 1;
    EXPECT_EQ(expected_result,  agent.bulk_contains(24));
    std::filesystem::remove(tmp_dir/"Test_IBF_1");
    std::filesystem::remove(tmp_dir/"Test_IBF_2");
}

TEST(ibf, no_given_expression_levels)
{
    arguments args{};
    ibf_arguments ibf_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    ibf_args.set_expression_levels_samplewise = true;
    ibf_args.number_expression_levels = 2;
    ibf_args.sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};

    std::vector<uint32_t> expected{2, 3};

    std::vector<uint32_t> medians = ibf(args, ibf_args);

    EXPECT_EQ(expected, medians);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    load_ibf(ibf, tmp_dir/"Test_IBF_Level_1");
    auto agent = ibf.membership_agent();

    sdsl::bit_vector expected_result(1, 0);
    EXPECT_EQ(expected_result,  agent.bulk_contains(2));
    expected_result[0] = 1;
    EXPECT_EQ(expected_result,  agent.bulk_contains(97));
    std::filesystem::remove(tmp_dir/"Test_IBF_Level_1");
}

TEST(ibf, no_given_expression_levels_auto)
{
    arguments args{};
    ibf_arguments ibf_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    ibf_args.number_expression_levels = 2;
    ibf_args.sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};

    std::vector<uint32_t> expected{2, 4};

    std::vector<uint32_t> medians = ibf(args, ibf_args);

    EXPECT_EQ(expected, medians);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    load_ibf(ibf, tmp_dir/"Test_IBF_2");
    auto agent = ibf.membership_agent();

    sdsl::bit_vector expected_result(1, 0);
    EXPECT_EQ(expected_result,  agent.bulk_contains(2));
    expected_result[0] = 1;
    EXPECT_EQ(expected_result,  agent.bulk_contains(97));
    std::filesystem::remove(tmp_dir/"Test_IBF_2");
    std::filesystem::remove(tmp_dir/"Test_IBF_4");
}

TEST(ibf, throws)
{
    arguments args{};
    ibf_arguments ibf_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    ibf_args.set_expression_levels_samplewise = true;
    ibf_args.sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};

    EXPECT_THROW(ibf(args, ibf_args), std::invalid_argument);

    ibf_args.expression_levels = {1, 2};
    ibf_args.bin_size = {};
    EXPECT_THROW(ibf(args, ibf_args), std::invalid_argument);

    ibf_args.bin_size = {1000};
    EXPECT_THROW(ibf(args, ibf_args), std::invalid_argument);
}

TEST(ibfmin, given_expression_levels)
{
    arguments args{};
    ibf_arguments ibf_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    ibf_args.expression_levels = {1, 2};
    ibf_args.bin_size = {1000, 1000};
    std::vector<std::filesystem::path> minimiser_file = {std::string(DATA_INPUT_DIR) + "mini_example.minimiser"};

    std::vector<uint32_t> expected{1, 2};

    std::vector<uint32_t> medians = ibf(minimiser_file, args, ibf_args);

    EXPECT_EQ(expected, medians);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    load_ibf(ibf, tmp_dir/"Test_IBF_1");
    auto agent = ibf.membership_agent();

    sdsl::bit_vector expected_result(1, 0);
    EXPECT_EQ(expected_result,  agent.bulk_contains(97));
    expected_result[0] = 1;
    EXPECT_EQ(expected_result,  agent.bulk_contains(24));
    std::filesystem::remove(tmp_dir/"Test_IBF_1");
    std::filesystem::remove(tmp_dir/"Test_IBF_2");
}


TEST(ibfmin, no_given_expression_levels)
{
    arguments args{};
    ibf_arguments ibf_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    ibf_args.set_expression_levels_samplewise = true;
    ibf_args.number_expression_levels = 2;
    ibf_args.bin_size = {1000, 1000};
    std::vector<std::filesystem::path> minimiser_file = {std::string(DATA_INPUT_DIR) + "mini_example.minimiser"};

    std::vector<uint32_t> expected{2, 3};

    std::vector<uint32_t> medians = ibf(minimiser_file, args, ibf_args);

    EXPECT_EQ(expected, medians);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    load_ibf(ibf, tmp_dir/"Test_IBF_Level_0");
    auto agent = ibf.membership_agent();

    sdsl::bit_vector expected_result(1, 0);
    EXPECT_EQ(expected_result,  agent.bulk_contains(2));
    expected_result[0] = 1;
    EXPECT_EQ(expected_result,  agent.bulk_contains(108));
    std::filesystem::remove(tmp_dir/"Test_IBF_Level_0");
    std::filesystem::remove(tmp_dir/"Test_IBF_Levels.levels");
}

TEST(ibfmin, no_given_expression_levels_auto)
{
    arguments args{};
    ibf_arguments ibf_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    ibf_args.number_expression_levels = 2;
    ibf_args.bin_size = {1000, 1000};
    std::vector<std::filesystem::path> minimiser_file = {std::string(DATA_INPUT_DIR) + "mini_example.minimiser"};

    std::vector<uint32_t> expected{2, 4};

    std::vector<uint32_t> medians = ibf(minimiser_file, args, ibf_args);

    EXPECT_EQ(expected, medians);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    load_ibf(ibf, tmp_dir/"Test_IBF_2");
    auto agent = ibf.membership_agent();

    sdsl::bit_vector expected_result(1, 0);
    EXPECT_EQ(expected_result,  agent.bulk_contains(2));
    expected_result[0] = 1;
    EXPECT_EQ(expected_result,  agent.bulk_contains(97));
    std::filesystem::remove(tmp_dir/"Test_IBF_2");
    std::filesystem::remove(tmp_dir/"Test_IBF_4");
}

TEST(minimiser, small_example)
{
    arguments args{};
    ibf_arguments ibf_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    ibf_args.expression_levels = {0};
    ibf_args.sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta",
                               std::string(DATA_INPUT_DIR) + "mini_example2.fasta"};

    minimiser(args, ibf_args);
    std::vector<uint64_t> counts{};
    uint32_t normalized_exp_value{};
    robin_hood::unordered_node_map<uint64_t, uint16_t> result_hash_table{};
    std::vector<std::filesystem::path> minimiser_files{};
    seqan3::shape expected_shape = seqan3::ungapped{args.k};

    for (int i = 0; i < ibf_args.sequence_files.size(); ++i)
    {
        // Test Header file
        ibf_args.expression_levels = {1, 2};
        read_header(args, ibf_args, std::string{ibf_args.path_out}  +
                    std::string{ibf_args.sequence_files[i].stem()} + ".header", counts);

        EXPECT_EQ(4, args.k);
        EXPECT_EQ(4, args.w_size.get());
        EXPECT_EQ(0, args.s.get());
        EXPECT_EQ(15, args.shape.to_ulong());
        EXPECT_EQ(0, ibf_args.cutoffs[0]);
        EXPECT_EQ(0, ibf_args.expression_levels[2]);
        EXPECT_EQ(3, ibf_args.expression_levels.size());
        EXPECT_EQ(12, counts[0]);
        EXPECT_EQ(1, counts.size());
        counts.clear();

        // Test binary file
        read_binary(result_hash_table, tmp_dir/("Test_" + std::string{ibf_args.sequence_files[i].stem()} + ".minimiser"));
        minimiser_files.push_back(tmp_dir/("Test_" + std::string{ibf_args.sequence_files[i].stem()} + ".minimiser"));
        EXPECT_EQ(expected_hash_tables[i], result_hash_table);

        result_hash_table.clear();
    }

    EXPECT_EQ(ibf_args.expression_levels, ibf(minimiser_files, args, ibf_args));

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    load_ibf(ibf, tmp_dir/"Test_IBF_0");
    auto agent = ibf.membership_agent();

    sdsl::bit_vector expected_result(2, 0);
    EXPECT_EQ(expected_result,  agent.bulk_contains(2));
    expected_result[0] = 1;
    EXPECT_EQ(expected_result,  agent.bulk_contains(0));
    expected_result[1] = 1;
    EXPECT_EQ(expected_result,  agent.bulk_contains(27));

    std::filesystem::remove(tmp_dir/"Test_IBF_0");
    std::filesystem::remove(tmp_dir/("Test_mini_example.header"));
    std::filesystem::remove(tmp_dir/("Test_mini_example2.header"));
    std::filesystem::remove(tmp_dir/("Test_mini_example.minimiser"));
    std::filesystem::remove(tmp_dir/("Test_mini_example2.minimiser"));
}

TEST(minimiser, small_example_auto_expression_level)
{
    arguments args{};
    ibf_arguments ibf_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    ibf_args.number_expression_levels = 2;
    ibf_args.sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta",
                               std::string(DATA_INPUT_DIR) + "mini_example2.fasta"};

    minimiser(args, ibf_args);
    std::vector<uint64_t> counts{};
    uint32_t normalized_exp_value{};
    std::vector<std::vector<uint64_t>> expected_counts{{6, 3}, {1, 1}};
    std::vector<uint64_t> expected_levels{2, 4};
    robin_hood::unordered_node_map<uint64_t, uint16_t> result_hash_table{};
    std::vector<std::filesystem::path> minimiser_files{};
    seqan3::shape expected_shape = seqan3::ungapped{args.k};

    for (int i = 0; i < ibf_args.sequence_files.size(); ++i)
    {
        // Test Header file
        ibf_args.expression_levels = {};
        read_header(args, ibf_args, std::string{ibf_args.path_out}  +
                    std::string{ibf_args.sequence_files[i].stem()} + ".header", counts);
        EXPECT_EQ(4, args.k);
        EXPECT_EQ(4, args.w_size.get());
        EXPECT_EQ(0, args.s.get());
        EXPECT_EQ(15, args.shape.to_ulong());
        EXPECT_EQ(0, ibf_args.cutoffs[0]);
        EXPECT_EQ(expected_levels[0], ibf_args.expression_levels[0]);
        EXPECT_EQ(expected_levels[1], ibf_args.expression_levels[1]);
        EXPECT_EQ(2, ibf_args.expression_levels.size());
        EXPECT_EQ(expected_counts[i][0], counts[0]);
        EXPECT_EQ(expected_counts[i][1], counts[1]);
        EXPECT_EQ(2, counts.size());
        counts.clear();

        // Test binary file
        read_binary(result_hash_table, tmp_dir/("Test_" + std::string{ibf_args.sequence_files[i].stem()} + ".minimiser"));
        minimiser_files.push_back(tmp_dir/("Test_" + std::string{ibf_args.sequence_files[i].stem()} + ".minimiser"));
        EXPECT_EQ(expected_hash_tables[i], result_hash_table);

        result_hash_table.clear();
    }

    EXPECT_EQ(ibf_args.expression_levels, ibf(minimiser_files, args, ibf_args));

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    load_ibf(ibf, tmp_dir/"Test_IBF_2");
    auto agent = ibf.membership_agent();

    sdsl::bit_vector expected_result(2, 0);
    EXPECT_EQ(expected_result, agent.bulk_contains(2));
    EXPECT_EQ(expected_result, agent.bulk_contains(27));
    expected_result[0] = 1;
    EXPECT_EQ(expected_result, agent.bulk_contains(0));
    std::filesystem::remove(tmp_dir/"Test_IBF_2");
    std::filesystem::remove(tmp_dir/("Test_mini_example.header"));
    std::filesystem::remove(tmp_dir/("Test_mini_example2.header"));
    std::filesystem::remove(tmp_dir/("Test_mini_example.minimiser"));
    std::filesystem::remove(tmp_dir/("Test_mini_example2.minimiser"));
}

TEST(minimiser, small_example_samplewise)
{
    arguments args{};
    ibf_arguments ibf_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);

    ibf_args.number_expression_levels = 1;
    ibf_args.sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta",
                               std::string(DATA_INPUT_DIR) + "mini_example2.fasta"};
    ibf_args.set_expression_levels_samplewise = true;

    minimiser(args, ibf_args);
    std::vector<uint64_t> counts{};
    uint32_t normalized_exp_value{};
    std::vector<std::vector<uint64_t>> expected_counts{{7}, {12}};
    std::vector<uint64_t> expected_levels{3, 1};
    robin_hood::unordered_node_map<uint64_t, uint16_t> result_hash_table{};
    std::vector<std::filesystem::path> minimiser_files{};
    seqan3::shape expected_shape = seqan3::ungapped{args.k};

    for (int i = 0; i < ibf_args.sequence_files.size(); ++i)
    {
        // Test Header file
        ibf_args.expression_levels = {};
        read_header(args, ibf_args, std::string{ibf_args.path_out}  +
                    std::string{ibf_args.sequence_files[i].stem()} + ".header", counts);
        EXPECT_EQ(4, args.k);
        EXPECT_EQ(4, args.w_size.get());
        EXPECT_EQ(0, args.s.get());
        EXPECT_EQ(15, args.shape.to_ulong());
        EXPECT_EQ(0, ibf_args.cutoffs[0]);
        EXPECT_EQ(expected_levels[i], ibf_args.expression_levels[0]);
        EXPECT_EQ(1, ibf_args.expression_levels.size());
        EXPECT_EQ(expected_counts[i][0], counts[0]);
        EXPECT_EQ(1, counts.size());
        counts.clear();

        // Test binary file
        read_binary(result_hash_table, tmp_dir/("Test_" + std::string{ibf_args.sequence_files[i].stem()} + ".minimiser"));
        minimiser_files.push_back(tmp_dir/("Test_" + std::string{ibf_args.sequence_files[i].stem()} + ".minimiser"));
        EXPECT_EQ(expected_hash_tables[i], result_hash_table);

        result_hash_table.clear();
        expected_levels = {1, 1};
    }

    EXPECT_EQ(ibf_args.expression_levels, ibf(minimiser_files, args, ibf_args));

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    load_ibf(ibf, tmp_dir/"Test_IBF_Level_0");
    auto agent = ibf.membership_agent();

    sdsl::bit_vector expected_result(2, 0);
    EXPECT_EQ(expected_result, agent.bulk_contains(2));
    EXPECT_EQ(expected_result, agent.bulk_contains(0));
    expected_result[0] = 1;
    expected_result[1] = 1;
    EXPECT_EQ(expected_result, agent.bulk_contains(27));
    std::filesystem::remove(tmp_dir/"Test_IBF_Level_0");
    std::filesystem::remove(tmp_dir/("Test_mini_example.header"));
    std::filesystem::remove(tmp_dir/("Test_mini_example2.header"));
    std::filesystem::remove(tmp_dir/("Test_mini_example.minimiser"));
    std::filesystem::remove(tmp_dir/("Test_mini_example2.minimiser"));
}

TEST(estimate, small_example)
{
    arguments args{};
    ibf_arguments ibf_args{};
    estimate_arguments estimate_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    ibf_args.expression_levels = {1, 2, 4};
    ibf_args.sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};
    estimate_args.expressions = ibf_args.expression_levels;

    ibf(args, ibf_args);
    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    estimate(args, estimate_args, ibf, tmp_dir/"expression.out",
             std::string(DATA_INPUT_DIR) + "mini_gen.fasta", ibf_args.path_out);

    std::ifstream output_file(tmp_dir/"expression.out");
    std::string line;
    std::string expected{"gen1\t3\t"};
    if (output_file.is_open())
    {
        while ( std::getline (output_file,line) )
        {
            EXPECT_EQ(expected,line);
        }
        output_file.close();
    }
    std::filesystem::remove(tmp_dir/"Test_IBF_2");
    std::filesystem::remove(tmp_dir/"Test_IBF_4");
    std::filesystem::remove(tmp_dir/"expression.out");
}

TEST(estimate, small_example_uncompressed)
{
    arguments args{};
    ibf_arguments ibf_args{};
    estimate_arguments estimate_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    args.compressed = false;
    ibf_args.expression_levels = {1, 2, 4};
    ibf_args.sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};
    estimate_args.expressions = ibf_args.expression_levels;

    ibf(args, ibf_args);
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf;
    estimate(args, estimate_args, ibf, tmp_dir/"expression.out",
             std::string(DATA_INPUT_DIR) + "mini_gen.fasta", ibf_args.path_out);

    std::ifstream output_file(tmp_dir/"expression.out");
    std::string line;
    std::string expected{"gen1\t3\t"};
    if (output_file.is_open())
    {
        while ( std::getline (output_file,line) )
        {
            EXPECT_EQ(expected,line);
        }
        output_file.close();
    }
    std::filesystem::remove(tmp_dir/"Test_IBF_2");
    std::filesystem::remove(tmp_dir/"Test_IBF_4");
    std::filesystem::remove(tmp_dir/"expression.out");
}

TEST(estimate, small_example_gene_not_found)
{
    arguments args{};
    ibf_arguments ibf_args{};
    estimate_arguments estimate_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    ibf_args.expression_levels = {2, 4};
    estimate_args.expressions = ibf_args.expression_levels;
    ibf_args.sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};

    ibf(args, ibf_args);
    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    estimate(args, estimate_args, ibf, tmp_dir/"expression.out",
             std::string(DATA_INPUT_DIR) + "mini_gen2.fasta", ibf_args.path_out);

    std::ifstream output_file(tmp_dir/"expression.out");
    std::string line;
    std::string expected{"gen2\t0\t"};
    if (output_file.is_open())
    {
        while ( std::getline (output_file,line) )
        {
            EXPECT_EQ(expected,line);
        }
        output_file.close();
    }
    std::filesystem::remove(tmp_dir/"Test_IBF_2");
    std::filesystem::remove(tmp_dir/"Test_IBF_4");
    std::filesystem::remove(tmp_dir/"expression.out");
}

TEST(estimate, small_example_different_expressions_per_level)
{
    arguments args{};
    ibf_arguments ibf_args{};
    estimate_arguments estimate_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    ibf_args.number_expression_levels = 3;
    ibf_args.set_expression_levels_samplewise = true;
    ibf_args.sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};

    minimiser(args, ibf_args);
    std::vector<std::filesystem::path> minimiser_files{std::string(DATA_INPUT_DIR) + "mini_example.minimiser"};
    ibf_args.expression_levels = {};
    ibf(minimiser_files, args, ibf_args);
    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    estimate_args.expressions = {0, 1, 2};
    estimate(args, estimate_args, ibf, tmp_dir/"expression.out",
             std::string(DATA_INPUT_DIR) + "mini_gen.fasta", ibf_args.path_out, tmp_dir/"Test_IBF_Levels.levels");

    std::ifstream output_file(tmp_dir/"expression.out");
    std::string line;
    std::string expected{"gen1\t3\t"};
    if (output_file.is_open())
    {
        while ( std::getline (output_file,line) )
        {
            EXPECT_EQ(expected,line);
        }
        output_file.close();
    }
    std::filesystem::remove(tmp_dir/"Test_IBF_Level_0");
    std::filesystem::remove(tmp_dir/"Test_IBF_Level_1");
    std::filesystem::remove(tmp_dir/"Test_IBF_Level_2");
    std::filesystem::remove(tmp_dir/"Test_IBF_Levels.levels");
    std::filesystem::remove(tmp_dir/"expression.out");
}

TEST(estimate, example)
{
    arguments args{};
    ibf_arguments ibf_args{};
    estimate_arguments estimate_args{};
    ibf_args.sequence_files = {std::string(DATA_INPUT_DIR) + "exp_01.fasta", std::string(DATA_INPUT_DIR) + "exp_02.fasta",
                               std::string(DATA_INPUT_DIR) + "exp_11.fasta", std::string(DATA_INPUT_DIR) + "exp_12.fasta"};
    ibf_args.samples = {2,2};
    ibf_args.expression_levels = {32};
    ibf_args.bin_size = {100000};
    ibf_args.path_out = tmp_dir/"";
    args.compressed = false;
    estimate_args.expressions = ibf_args.expression_levels;
    ibf(args, ibf_args);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf;
    estimate(args, estimate_args, ibf, tmp_dir/"expression.out",
    std::string(DATA_INPUT_DIR) + "gene.fasta", ibf_args.path_out);

    std::ifstream output_file(tmp_dir/"expression.out");
    std::string line;
    std::string expected{"GeneA\t0\t32\t"};
    if (output_file.is_open())
    {
        while ( std::getline (output_file,line) )
        {
             EXPECT_EQ(expected,line);
        }
        output_file.close();
    }
    std::filesystem::remove(tmp_dir/"Test_IBF_32");
    std::filesystem::remove(tmp_dir/"expression.out");
}

TEST(estimate, example_different_expressions_per_level)
{
    arguments args{};
    ibf_arguments ibf_args{};
    estimate_arguments estimate_args{};
    ibf_args.cutoffs = {1};
    ibf_args.sequence_files = {std::string(DATA_INPUT_DIR) + "exp_01.fasta", std::string(DATA_INPUT_DIR) + "exp_02.fasta",
                               std::string(DATA_INPUT_DIR) + "exp_11.fasta", std::string(DATA_INPUT_DIR) + "exp_12.fasta"};
    ibf_args.samples = {2,2};
    ibf_args.number_expression_levels = 3;
    ibf_args.set_expression_levels_samplewise = true;
    ibf_args.bin_size = {10000000};
    ibf_args.path_out = tmp_dir/"Test_";
    args.compressed = false;
    estimate_args.expressions = {0, 1, 2};
    minimiser(args, ibf_args);
    std::vector<std::filesystem::path> minimiser_files{std::string(DATA_INPUT_DIR) + "exp_01.minimiser", std::string(DATA_INPUT_DIR) + "exp_11.minimiser"};
    ibf_args.expression_levels = {};
    ibf(minimiser_files, args, ibf_args);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf;
    estimate(args, estimate_args, ibf, tmp_dir/"expression.out",
    std::string(DATA_INPUT_DIR) + "gene.fasta", ibf_args.path_out, tmp_dir/"Test_IBF_Levels.levels");

    std::ifstream output_file(tmp_dir/"expression.out");
    std::string line;
    // Count would expect 6 and 44
    std::string expected{"GeneA\t8\t29\t"};
    if (output_file.is_open())
    {
        while ( std::getline (output_file,line) )
        {
             EXPECT_EQ(expected,line);
        }
        output_file.close();
    }
    std::filesystem::remove(tmp_dir/"Test_IBF_Level_0");
    std::filesystem::remove(tmp_dir/"Test_IBF_Level_1");
    std::filesystem::remove(tmp_dir/"Test_IBF_Level_2");
    std::filesystem::remove(tmp_dir/"Test_IBF_Levels.levels");
    std::filesystem::remove(tmp_dir/"expression.out");
}

TEST(stats, example)
{
    arguments args{};
    ibf_arguments ibf_args{};
    ibf_args.cutoffs = {1};
    std::vector<std::filesystem::path> minimiser_files{std::string(DATA_INPUT_DIR) + "exp_01.header",
                                                       std::string(DATA_INPUT_DIR) + "exp_11.header"};

    std::vector<std::tuple<std::vector<uint64_t>, std::vector<uint64_t>>> expected{{{2}, {5310, 5739, 5739}},
                                                                                {{21}, {5235, 5904, 5904}},
                                                                                {{29}, {6115, 6358, 6358}}};

    std::vector<std::tuple<std::vector<uint64_t>, std::vector<uint64_t>>> results = statistics(args, ibf_args, minimiser_files);

    EXPECT_EQ(expected, results);
}
