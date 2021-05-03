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
    args.path_out = tmp_dir/"Test_";
}

void initialization_ibf_args(ibf_arguments & args)
{
    args.bin_size = {1000};
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
    minimiser_arguments minimiser_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    ibf_args.expression_levels = {1, 2};
    minimiser_args.sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};

    std::vector<uint16_t> expected{1, 2};

    std::vector<uint16_t> medians = ibf(args, ibf_args, minimiser_args);

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
    minimiser_arguments minimiser_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    ibf_args.number_expression_levels = 2;
    minimiser_args.sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};

    std::vector<uint16_t> expected{3, 4};

    std::vector<uint16_t> medians = ibf(args, ibf_args, minimiser_args);

    EXPECT_EQ(expected, medians);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    load_ibf(ibf, tmp_dir/"Test_IBF_Level_0");
    auto agent = ibf.membership_agent();

    sdsl::bit_vector expected_result(1, 0);
    EXPECT_EQ(expected_result,  agent.bulk_contains(2));
    expected_result[0] = 1;
    EXPECT_EQ(expected_result,  agent.bulk_contains(97));
    std::filesystem::remove(tmp_dir/"Test_IBF_Level_0");
}

TEST(ibf, throws)
{
    arguments args{};
    ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    minimiser_args.sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};

    EXPECT_THROW(ibf(args, ibf_args, minimiser_args), std::invalid_argument);

    ibf_args.number_expression_levels = 0;
    ibf_args.bin_size = {};
    EXPECT_THROW(ibf(args, ibf_args, minimiser_args), std::invalid_argument);

    ibf_args.bin_size = {1000};
    EXPECT_THROW(ibf(args, ibf_args, minimiser_args), std::invalid_argument);
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

    std::vector<uint16_t> expected{1, 2};

    std::vector<uint16_t> medians = ibf(minimiser_file, args, ibf_args);

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

TEST(ibfmin, given_expression_levels_multiple_threads)
{
    arguments args{};
    ibf_arguments ibf_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    ibf_args.expression_levels = {1, 2};
    ibf_args.bin_size = {1000, 1000};
    args.threads = 2;
    std::vector<std::filesystem::path> minimiser_file = {std::string(DATA_INPUT_DIR) + "mini_example.minimiser", std::string(DATA_INPUT_DIR) + "mini_example.minimiser"};

    std::vector<uint16_t> expected{1, 2};

    std::vector<uint16_t> medians = ibf(minimiser_file, args, ibf_args);

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
    ibf_args.number_expression_levels = 2;
    ibf_args.bin_size = {1000, 1000};
    std::vector<std::filesystem::path> minimiser_file = {std::string(DATA_INPUT_DIR) + "mini_example.minimiser"};

    std::vector<uint16_t> expected{};

    std::vector<uint16_t> medians = ibf(minimiser_file, args, ibf_args);

    EXPECT_EQ(expected, medians);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    load_ibf(ibf, tmp_dir/"Test_IBF_Level_0");
    auto agent = ibf.membership_agent();

    sdsl::bit_vector expected_result(1, 0);
    EXPECT_EQ(expected_result,  agent.bulk_contains(2));
    expected_result[0] = 1;
    EXPECT_EQ(expected_result,  agent.bulk_contains(97));
    std::filesystem::remove(tmp_dir/"Test_IBF_Level_0");
    std::filesystem::remove(tmp_dir/"Test_IBF_Levels.levels");
}

TEST(ibfmin, no_given_expression_levels_multiple_threads)
{
    arguments args{};
    ibf_arguments ibf_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    ibf_args.number_expression_levels = 2;
    ibf_args.bin_size = {1000, 1000};
    args.threads = 2;
    std::vector<std::filesystem::path> minimiser_file = {std::string(DATA_INPUT_DIR) + "mini_example.minimiser", std::string(DATA_INPUT_DIR) + "mini_example.minimiser"};

    std::vector<uint16_t> expected{};

    std::vector<uint16_t> medians = ibf(minimiser_file, args, ibf_args);

    EXPECT_EQ(expected, medians);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    load_ibf(ibf, tmp_dir/"Test_IBF_Level_0");
    auto agent = ibf.membership_agent();

    sdsl::bit_vector expected_result(1, 0);
    EXPECT_EQ(expected_result,  agent.bulk_contains(2));
    expected_result[0] = 1;
    EXPECT_EQ(expected_result,  agent.bulk_contains(97));
    std::filesystem::remove(tmp_dir/"Test_IBF_Level_0");
    std::filesystem::remove(tmp_dir/"Test_IBF_Levels.levels");
}

TEST(minimiser, small_example)
{
    arguments args{};
    ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    minimiser_args.cutoffs = {0, 0};
    ibf_args.expression_levels = {0};
    minimiser_args.sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta",
                               std::string(DATA_INPUT_DIR) + "mini_example2.fasta"};

    minimiser(args, minimiser_args);
    uint32_t normalized_exp_value{};
    robin_hood::unordered_node_map<uint64_t, uint16_t> result_hash_table{};
    std::vector<std::filesystem::path> minimiser_files{};
    seqan3::shape expected_shape = seqan3::ungapped{args.k};

    for (int i = 0; i < minimiser_args.sequence_files.size(); ++i)
    {
        // Test Header file
        read_header(args, minimiser_args.cutoffs, std::string{args.path_out}  +
                    std::string{minimiser_args.sequence_files[i].stem()} + ".header");

        EXPECT_EQ(4, args.k);
        EXPECT_EQ(4, args.w_size.get());
        EXPECT_EQ(0, args.s.get());
        EXPECT_EQ(15, args.shape.to_ulong());
        EXPECT_EQ(0, minimiser_args.cutoffs[0]);

        // Test binary file
        read_binary(result_hash_table, tmp_dir/("Test_" + std::string{minimiser_args.sequence_files[i].stem()} + ".minimiser"));
        minimiser_files.push_back(tmp_dir/("Test_" + std::string{minimiser_args.sequence_files[i].stem()} + ".minimiser"));
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

TEST(minimiser, small_example_samplewise)
{
    arguments args{};
    ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);

    minimiser_args.cutoffs = {0, 0};
    ibf_args.number_expression_levels = 1;
    minimiser_args.sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta",
                               std::string(DATA_INPUT_DIR) + "mini_example2.fasta"};

    minimiser(args, minimiser_args);
    uint32_t normalized_exp_value{};
    std::vector<std::vector<uint32_t>> expected_counts{{7}, {12}};
    std::vector<uint16_t> expected_levels{3, 1};
    robin_hood::unordered_node_map<uint64_t, uint16_t> result_hash_table{};
    std::vector<std::filesystem::path> minimiser_files{};
    seqan3::shape expected_shape = seqan3::ungapped{args.k};

    for (int i = 0; i < minimiser_args.sequence_files.size(); ++i)
    {
        // Test Header file
        ibf_args.expression_levels = {};
        read_header(args, minimiser_args.cutoffs, std::string{args.path_out}  +
                    std::string{minimiser_args.sequence_files[i].stem()} + ".header");
        EXPECT_EQ(4, args.k);
        EXPECT_EQ(4, args.w_size.get());
        EXPECT_EQ(0, args.s.get());
        EXPECT_EQ(15, args.shape.to_ulong());
        EXPECT_EQ(0, minimiser_args.cutoffs[0]);

        // Test binary file
        read_binary(result_hash_table, tmp_dir/("Test_" + std::string{minimiser_args.sequence_files[i].stem()} + ".minimiser"));
        minimiser_files.push_back(tmp_dir/("Test_" + std::string{minimiser_args.sequence_files[i].stem()} + ".minimiser"));
        EXPECT_EQ(expected_hash_tables[i], result_hash_table);

        result_hash_table.clear();
    }
    ibf_args.expression_levels = {};
    expected_levels = {}; // Only levels from last experiment are returned.
    EXPECT_EQ(expected_levels, ibf(minimiser_files, args, ibf_args));

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

TEST(minimiser, cutoff_by_filesize)
{
    arguments args{};
    ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    ibf_args.expression_levels = {0};
    minimiser_args.sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta",
                               std::string(DATA_INPUT_DIR) + "mini_example2.fasta"};

    minimiser(args, minimiser_args);
    uint32_t normalized_exp_value{};
    std::vector<std::filesystem::path> minimiser_files{};
    seqan3::shape expected_shape = seqan3::ungapped{args.k};

    for (int i = 0; i < minimiser_args.sequence_files.size(); ++i)
    {
        // Test Header file
        read_header(args, minimiser_args.cutoffs, std::string{args.path_out}  +
                    std::string{minimiser_args.sequence_files[i].stem()} + ".header");

        EXPECT_EQ(4, args.k);
        EXPECT_EQ(4, args.w_size.get());
        EXPECT_EQ(0, args.s.get());
        EXPECT_EQ(15, args.shape.to_ulong());
        EXPECT_EQ(1, minimiser_args.cutoffs[0]);
        minimiser_files.push_back(tmp_dir/("Test_" + std::string{minimiser_args.sequence_files[i].stem()} + ".minimiser"));
    }

    EXPECT_EQ(ibf_args.expression_levels, ibf(minimiser_files, args, ibf_args));

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    load_ibf(ibf, tmp_dir/"Test_IBF_0");
    auto agent = ibf.membership_agent();

    sdsl::bit_vector expected_result(2, 0);
    EXPECT_EQ(expected_result,  agent.bulk_contains(2));
    expected_result[0] = 1;
    EXPECT_EQ(expected_result,  agent.bulk_contains(0));
    expected_result[0] = 0;
    expected_result[1] = 1;
    EXPECT_EQ(expected_result,  agent.bulk_contains(85));

    std::filesystem::remove(tmp_dir/"Test_IBF_0");
    std::filesystem::remove(tmp_dir/("Test_mini_example.header"));
    std::filesystem::remove(tmp_dir/("Test_mini_example2.header"));
    std::filesystem::remove(tmp_dir/("Test_mini_example.minimiser"));
    std::filesystem::remove(tmp_dir/("Test_mini_example2.minimiser"));
}

TEST(minimiser, small_example_two_threads)
{
    arguments args{};
    ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    args.threads = 2;
    minimiser_args.cutoffs = {0, 0};
    ibf_args.expression_levels = {0};
    minimiser_args.sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta",
                               std::string(DATA_INPUT_DIR) + "mini_example2.fasta"};
    minimiser(args, minimiser_args);
    args.threads = 1;
    uint32_t normalized_exp_value{};
    robin_hood::unordered_node_map<uint64_t, uint16_t> result_hash_table{};
    std::vector<std::filesystem::path> minimiser_files{};
    seqan3::shape expected_shape = seqan3::ungapped{args.k};

    for (int i = 0; i < minimiser_args.sequence_files.size(); ++i)
    {
        // Test Header file
        read_header(args, minimiser_args.cutoffs, std::string{args.path_out}  +
                    std::string{minimiser_args.sequence_files[i].stem()} + ".header");

        EXPECT_EQ(4, args.k);
        EXPECT_EQ(4, args.w_size.get());
        EXPECT_EQ(0, args.s.get());
        EXPECT_EQ(15, args.shape.to_ulong());
        EXPECT_EQ(0, minimiser_args.cutoffs[0]);

        // Test binary file
        read_binary(result_hash_table, tmp_dir/("Test_" + std::string{minimiser_args.sequence_files[i].stem()} + ".minimiser"));
        minimiser_files.push_back(tmp_dir/("Test_" + std::string{minimiser_args.sequence_files[i].stem()} + ".minimiser"));
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

TEST(estimate, small_example)
{
    arguments args{};
    ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    estimate_arguments estimate_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    ibf_args.expression_levels = {1, 2, 4};
    minimiser_args.sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};
    estimate_args.expressions = ibf_args.expression_levels;

    ibf(args, ibf_args, minimiser_args);
    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    estimate(args, estimate_args, ibf, tmp_dir/"expression.out",
             std::string(DATA_INPUT_DIR) + "mini_gen.fasta", args.path_out);

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
    minimiser_arguments minimiser_args{};
    estimate_arguments estimate_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    args.compressed = false;
    ibf_args.expression_levels = {1, 2, 4};
    minimiser_args.sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};
    estimate_args.expressions = ibf_args.expression_levels;

    ibf(args, ibf_args, minimiser_args);
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf;
    estimate(args, estimate_args, ibf, tmp_dir/"expression.out",
             std::string(DATA_INPUT_DIR) + "mini_gen.fasta", args.path_out);

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
    minimiser_arguments minimiser_args{};
    estimate_arguments estimate_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    ibf_args.expression_levels = {2, 4};
    estimate_args.expressions = ibf_args.expression_levels;
    minimiser_args.sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};

    ibf(args, ibf_args, minimiser_args);
    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    estimate(args, estimate_args, ibf, tmp_dir/"expression.out",
             std::string(DATA_INPUT_DIR) + "mini_gen2.fasta", args.path_out);

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
    minimiser_arguments minimiser_args{};
    estimate_arguments estimate_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    ibf_args.number_expression_levels = 3;
    minimiser_args.sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};

    minimiser(args, minimiser_args);
    std::vector<std::filesystem::path> minimiser_files{tmp_dir/"Test_mini_example.minimiser"};
    ibf_args.expression_levels = {};
    ibf(minimiser_files, args, ibf_args);
    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    estimate_args.expressions = {0, 1, 2};
    estimate(args, estimate_args, ibf, tmp_dir/"expression.out",
             std::string(DATA_INPUT_DIR) + "mini_gen.fasta", args.path_out, tmp_dir/"Test_IBF_Levels.levels");

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
    minimiser_arguments minimiser_args{};
    estimate_arguments estimate_args{};
    minimiser_args.sequence_files = {std::string(DATA_INPUT_DIR) + "exp_01.fasta", std::string(DATA_INPUT_DIR) + "exp_02.fasta",
                               std::string(DATA_INPUT_DIR) + "exp_11.fasta", std::string(DATA_INPUT_DIR) + "exp_12.fasta"};
    minimiser_args.samples = {2,2};
    ibf_args.expression_levels = {32};
    ibf_args.bin_size = {100000};
    args.path_out = tmp_dir/"Test_";
    args.compressed = false;
    estimate_args.expressions = ibf_args.expression_levels;
    ibf(args, ibf_args, minimiser_args);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf;
    estimate(args, estimate_args, ibf, tmp_dir/"expression.out",
    std::string(DATA_INPUT_DIR) + "gene.fasta", args.path_out);

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
    minimiser_arguments minimiser_args{};
    estimate_arguments estimate_args{};
    minimiser_args.sequence_files = {std::string(DATA_INPUT_DIR) + "exp_01.fasta", std::string(DATA_INPUT_DIR) + "exp_02.fasta",
                               std::string(DATA_INPUT_DIR) + "exp_11.fasta", std::string(DATA_INPUT_DIR) + "exp_12.fasta"};
    minimiser_args.cutoffs = {0, 0};
    minimiser_args.samples = {2,2};
    ibf_args.number_expression_levels = 3;
    ibf_args.bin_size = {10000000};
    args.path_out = tmp_dir/"Test_";
    args.compressed = false;
    estimate_args.expressions = {0, 1, 2};
    minimiser(args, minimiser_args);
    std::vector<std::filesystem::path> minimiser_files{tmp_dir/"Test_exp_01.minimiser", tmp_dir/"Test_exp_11.minimiser"};
    ibf_args.expression_levels = {};
    ibf(minimiser_files, args, ibf_args);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf;
    estimate(args, estimate_args, ibf, tmp_dir/"expression.out",
    std::string(DATA_INPUT_DIR) + "gene.fasta", args.path_out, tmp_dir/"Test_IBF_Levels.levels");

    std::ifstream output_file(tmp_dir/"expression.out");
    std::string line;
    // Count would expect 6 and 34
    std::string expected{"GeneA\t4\t26\t"};
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
