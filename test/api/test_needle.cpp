#include <gtest/gtest.h>
#include <iostream>

#include <robin_hood.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>

#include "ibf.h"
#include "minimiser.h"
#include "search.h"

#ifndef DATA_INPUT_DIR
#  define DATA_INPUT_DIR @DATA_INPUT_DIR@
#endif

using seqan3::operator""_shape;

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
    args.expression_levels = {1};
    args.bin_size = {1000};
    args.path_out = DATA_INPUT_DIR;
}

TEST(count, small_example)
{
    arguments args{};
    initialization_args(args);

    count(args, {std::string(DATA_INPUT_DIR) + "mini_example.fasta"}, std::string(DATA_INPUT_DIR) + "mini_gen.fasta",
          std::string(DATA_INPUT_DIR), false);

    std::ifstream output_file(std::string(DATA_INPUT_DIR) + "mini_example.count.out");
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
}

TEST(count, small_example_paired)
{
    arguments args{};
    initialization_args(args);

    count(args, {std::string(DATA_INPUT_DIR) + "mini_example.fasta", std::string(DATA_INPUT_DIR) + "mini_example.fasta"},
          std::string(DATA_INPUT_DIR) + "mini_gen.fasta",
          std::string(DATA_INPUT_DIR), true);

    std::ifstream output_file(std::string(DATA_INPUT_DIR) + "mini_example.count.out");
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
}

TEST(ibf, median)
{
    arguments args{};
    ibf_arguments ibf_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    ibf_args.sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};

    std::vector<uint32_t> expected{3};

    std::vector<uint32_t> medians = ibf(args, ibf_args);

    EXPECT_EQ(expected, medians);
}

TEST(ibfmin, median)
{
    arguments args{};
    ibf_arguments ibf_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    std::vector<std::filesystem::path> minimiser_file = {std::string(DATA_INPUT_DIR) + "mini_example.minimiser"};
    std::filesystem::path header_file = "";

    std::vector<uint32_t> expected{3};

    std::vector<uint32_t> medians = ibf(minimiser_file, header_file, args, ibf_args);

    EXPECT_EQ(expected, medians);
}

TEST(ibf, mean)
{
    arguments args{};
    ibf_arguments ibf_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    ibf_args.sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};
    ibf_args.normalization_method = "mean";

    std::vector<uint32_t> expected{2};

    std::vector<uint32_t> means = ibf(args, ibf_args);

    EXPECT_EQ(expected, means);
}

TEST(ibf, genom_median)
{
    arguments args{};
    ibf_arguments ibf_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    ibf_args.sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};
    ibf_args.include_file = std::string(DATA_INPUT_DIR) + "mini_genom.fasta";

    std::vector<uint32_t> expected{4};

    std::vector<uint32_t> medians = ibf(args, ibf_args);

    EXPECT_EQ(expected, medians);
}

// Test a genome which has so different minimisers that some minimisers in seq can not be found
TEST(ibf, genom_median_no_match)
{
    arguments args{};
    ibf_arguments ibf_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    ibf_args.sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example2.fasta"};
    ibf_args.include_file = std::string(DATA_INPUT_DIR) + "mini_genom.fasta";

    std::vector<uint32_t> expected{1};

    std::vector<uint32_t> medians = ibf(args, ibf_args);

    EXPECT_EQ(expected, medians);
}

TEST(ibf, genom_mean)
{
    arguments args{};
    ibf_arguments ibf_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    ibf_args.sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};
    ibf_args.include_file = std::string(DATA_INPUT_DIR) + "mini_genom.fasta";
    ibf_args.normalization_method = "mean";

    std::vector<uint32_t> expected{3};

    std::vector<uint32_t> means = ibf(args, ibf_args);

    EXPECT_EQ(expected, means);
}

TEST(minimiser, small_example)
{
    arguments args{};
    ibf_arguments ibf_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    std::vector<robin_hood::unordered_node_map<uint64_t,uint64_t>> expected_hash_tables{   // minimisers:
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
    std::vector<uint32_t> expected_normalized_exp_values{3,1};

    ibf_args.expression_levels = {0};
    ibf_args.sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta",
                               std::string(DATA_INPUT_DIR) + "mini_example2.fasta"};

    minimiser(args, ibf_args);
    std::vector<uint64_t> counts{};
    uint32_t normalized_exp_value{};
    robin_hood::unordered_node_map<uint64_t,uint64_t> result_hash_table{};
    std::vector<std::filesystem::path> minimiser_files{};
    seqan3::shape expected_shape = seqan3::ungapped{args.k};

    for (int i = 0; i < ibf_args.sequence_files.size(); ++i)
    {
        // Test Header file
        ibf_args.expression_levels = {};
        read_header(args, ibf_args, std::string{ibf_args.path_out}  +
                    std::string{ibf_args.sequence_files[i].stem()} + ".header", counts, normalized_exp_value);

        EXPECT_EQ(4, args.k);
        EXPECT_EQ(4, args.w_size.get());
        EXPECT_EQ(0, args.s.get());
        EXPECT_EQ(15, args.shape.to_ulong());
        EXPECT_EQ(0, ibf_args.cutoffs[0]);
        EXPECT_EQ(expected_normalized_exp_values[i], normalized_exp_value);
        EXPECT_EQ("median", ibf_args.normalization_method);
        EXPECT_EQ(0, ibf_args.expression_levels[0]);
        EXPECT_EQ(1, ibf_args.expression_levels.size());
        EXPECT_EQ(12, counts[0]);
        EXPECT_EQ(1, counts.size());
        counts.clear();

        // Test binary file
        read_binary(result_hash_table, std::string{ibf_args.path_out} + std::string{ibf_args.sequence_files[i].stem()}
                    + ".minimiser");
        minimiser_files.push_back(std::string{ibf_args.path_out} + std::string{ibf_args.sequence_files[i].stem()}
                                  + ".minimiser");
        EXPECT_EQ(expected_hash_tables[i], result_hash_table);

        result_hash_table.clear();
    }

    std::vector<uint32_t> medians = ibf(minimiser_files, "", args, ibf_args);
    EXPECT_EQ(expected_normalized_exp_values, medians);
}


TEST(estimate, small_example)
{
    arguments args{};
    ibf_arguments ibf_args{};
    search_arguments search_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    search_args.threshold = 0.5;
    ibf_args.expression_levels = {1, 2};
    ibf_args.sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};

    ibf(args, ibf_args);
    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf;
    estimate(args, search_args, ibf, ibf_args.expression_levels, std::string(DATA_INPUT_DIR) + "expression.out",
             std::string(DATA_INPUT_DIR) + "mini_gen.fasta", ibf_args.path_out);

    std::ifstream output_file(std::string(DATA_INPUT_DIR) + "expression.out");
    std::string line;
    std::string expected{"gen1\t2\t"};
    if (output_file.is_open())
    {
        while ( std::getline (output_file,line) )
        {
            EXPECT_EQ(expected,line);
        }
    output_file.close();
    }
}

TEST(search, small_example)
{
    arguments args{};
    ibf_arguments ibf_args{};
    search_arguments search_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    std::vector<uint32_t> expected{1};
    ibf_args.sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};

    ibf(args, ibf_args);

    search_args.search_file = std::string(DATA_INPUT_DIR) + "mini_gen.fasta";
    search_args.path_in = ibf_args.path_out;
    search_args.expression = 1;

    std::vector<uint32_t> results{search(args, search_args)};

    EXPECT_EQ(expected, results);
}

TEST(search, small_example_uncompressed)
{
    arguments args{};
    ibf_arguments ibf_args{};
    search_arguments search_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    std::vector<uint32_t> expected{1};
    ibf_args.sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};
    args.compressed = false;
    ibf(args, ibf_args);

    search_args.search_file = std::string(DATA_INPUT_DIR) + "mini_gen.fasta";
    search_args.path_in = ibf_args.path_out;
    search_args.expression = 1;

    std::vector<uint32_t> results{search(args, search_args)};

    EXPECT_EQ(expected, results);
}

TEST(search, small_example_gene_not_found)
{
    arguments args{};
    ibf_arguments ibf_args{};
    search_arguments search_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    std::vector<uint32_t> expected{0};
    ibf_args.sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};

    ibf(args, ibf_args);

    search_args.search_file = std::string(DATA_INPUT_DIR) + "mini_gen2.fasta";
    search_args.path_in = ibf_args.path_out;
    search_args.expression = 1;

    std::vector<uint32_t> results{search(args, search_args)};

    EXPECT_EQ(expected, results);
}

TEST(search, small_example_own_cutoffs)
{
    arguments args{};
    ibf_arguments ibf_args{};
    search_arguments search_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    std::vector<uint32_t> expected{0};
    ibf_args.sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};
    ibf_args.expression_levels = {0};
    ibf_args.cutoffs = {2};

    ibf(args, ibf_args);

    search_args.search_file = std::string(DATA_INPUT_DIR) + "mini_gen3.fasta";
    search_args.path_in = ibf_args.path_out;
    search_args.expression = 1;

    std::vector<uint32_t> results{search(args, search_args)};

    EXPECT_EQ(expected, results);
}

TEST(search, threshold)
{
    arguments args{};
    ibf_arguments ibf_args{};
    search_arguments search_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    std::vector<uint32_t> expected_05{1};
    std::vector<uint32_t> expected_07{0};
    ibf_args.sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};

    ibf(args, ibf_args);

    search_args.search_file = std::string(DATA_INPUT_DIR) + "mini_gen4.fasta";
    search_args.path_in = ibf_args.path_out;
    search_args.expression = 1;

    std::vector<uint32_t> results_05{search(args, search_args)};
    EXPECT_EQ(expected_05, results_05);

    search_args.threshold = 0.7;
    std::vector<uint32_t> results_07{search(args, search_args)};

    EXPECT_EQ(expected_07, results_07);
}

TEST(search, expression_file)
{
    arguments args{};
    ibf_arguments ibf_args{};
    search_arguments search_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    std::vector<uint32_t> expected{1};
    ibf_args.sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};
    ibf_args.expression_levels = {1};
    ibf_args.cutoffs = {2};

    ibf(args, ibf_args);

    search_args.search_file = std::string(DATA_INPUT_DIR) + "mini_genes.fasta";
    search_args.path_in = ibf_args.path_out;
    search_args.exp_file  = {std::string(DATA_INPUT_DIR) + "mini_example_exp.tsv"};

    std::vector<uint32_t> results{search(args, search_args)};

    EXPECT_EQ(expected, results);
}

TEST(search, example)
{
    arguments args{};
    ibf_arguments ibf_args{};
    search_arguments search_args{};
    std::vector<uint32_t> expected{0,1};
    ibf_args.sequence_files = {std::string(DATA_INPUT_DIR) + "exp_01.fasta", std::string(DATA_INPUT_DIR) + "exp_02.fasta",
                               std::string(DATA_INPUT_DIR) + "exp_11.fasta", std::string(DATA_INPUT_DIR) + "exp_12.fasta"};
    ibf_args.samples = {2,2};
    ibf_args.expression_levels = {32};
    ibf_args.bin_size = {100000};
    ibf_args.path_out = std::string(DATA_INPUT_DIR);
    args.compressed = false;
    ibf(args, ibf_args);

    // ./bin/needle search DATA_INPUT_DIR"+"/gene.fasta -i DATA_INPUT_DIR"+"/ -e 32
    search_args.search_file = std::string(DATA_INPUT_DIR) + "gene.fasta";
    search_args.path_in = ibf_args.path_out;
    search_args.expression = ibf_args.expression_levels[0];

    std::vector<uint32_t> results{search(args, search_args)};

    EXPECT_EQ(expected, results);
}

TEST(stats, example)
{
    arguments args{};
    ibf_arguments ibf_args{};
    std::vector<std::filesystem::path> minimiser_files{std::string(DATA_INPUT_DIR) + "exp_01.header",
                                                       std::string(DATA_INPUT_DIR) + "exp_11.header"};

    std::vector<std::tuple<std::vector<uint64_t>, std::vector<uint64_t>>> expected{{{0, 29}, {62496, 63053, 63053}},
                                                                                {{1, 29}, {6116, 6359, 6359}},
                                                                                {{4, 29}, {7, 25, 25}}};

    std::vector<std::tuple<std::vector<uint64_t>, std::vector<uint64_t>>> results = statistics(args, ibf_args, minimiser_files);

    EXPECT_EQ(expected, results);
}
