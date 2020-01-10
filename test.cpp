#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>

#include "ibf.h"
#include "minimizer3.h"
#include "search.h"

TEST(minimizer, small_example)
{
    using seqan3::operator""_dna4;

    std::vector<seqan3::dna4_vector> seq = {"ACGTCGACGTTTAG"_dna4};
    std::vector<uint32_t> expected{2,2,2,1,2,2,2,1,1,1,1};

    auto hash_table = compute_occurrences(seq, 4, 4, 0, 0);

    std::vector<uint32_t> minimizer_occurences;
    for (auto & minHash : compute_minimizer(seq[0], 4, 4, 0, 0))
    {
        minimizer_occurences.push_back(hash_table[minHash]);
    }

    EXPECT_EQ(expected,  minimizer_occurences);
}

TEST(minimizer, small_example_k_w_unequal)
{
    using seqan3::operator""_dna4;

    std::vector<seqan3::dna4_vector> seq = {"ACGTCGACGTTTAG"_dna4};
    std::vector<uint32_t> expected{2,1,2,1,1};

    auto hash_table = compute_occurrences(seq, 4, 8, 0, 0);

    std::vector<uint32_t> minimizer_occurences;
    for (auto & minHash : compute_minimizer(seq[0], 4, 8, 0, 0))
    {
        minimizer_occurences.push_back(hash_table[minHash]);
    }

    EXPECT_EQ(expected, minimizer_occurences);
}

TEST(minimizer, small_example_gaps)
{
    using seqan3::operator""_dna4;

    std::vector<seqan3::dna4_vector> seq = {"ACGTCGACGTTTAG"_dna4};
    std::vector<uint32_t> expected{2,1,2,1,1};

    auto hash_table = compute_occurrences(seq, 5, 8, 0b10101, 0);

    std::vector<uint32_t> minimizer_occurences;
    for (auto & minHash : compute_minimizer(seq[0], 5, 8, 0b10101, 0))
    {
        minimizer_occurences.push_back(hash_table[minHash]);
    }

    EXPECT_EQ(expected, minimizer_occurences);
}

TEST(ibf, median)
{
    cmd_arguments args{};
    args.sequence_files = {"./example/mini_example.fasta"};
    args.expression_levels = {1};
    args.bits = {1000};
    args.path_out = "./example/";
    args.compressed = true;
    args.k = 4;
    args.window_size = 4;
    args.seed = 0;

    std::vector<uint32_t> expected{3};

    std::vector<uint32_t> medians = ibf(args);

    EXPECT_EQ(expected, medians);
}

TEST(ibf, mean)
{
    cmd_arguments args{};
    args.sequence_files = {"./example/mini_example.fasta"};
    args.expression_levels = {1};
    args.bits = {1000};
    args.path_out = "./example/";
    args.aggregate_by = "mean";
    args.compressed = true;
    args.k = 4;
    args.window_size = 4;
    args.seed = 0;

    std::vector<uint32_t> expected{3};

    std::vector<uint32_t> means = ibf(args);

    EXPECT_EQ(expected, means);
}

TEST(ibf, random)
{
    cmd_arguments args{};
    args.sequence_files = {"./example/mini_example.fasta"};
    args.expression_levels = {1};
    args.bits = {1000};
    args.path_out = "./example/";
    args.aggregate_by = "random";
    args.random = 40;
    args.compressed = true;
    args.k = 4;
    args.window_size = 4;
    args.seed = 0;

    std::vector<uint32_t> medians = ibf(args);

    EXPECT_TRUE((3 >= medians[0]) & (2 <= medians[0]));
}

TEST(ibf, genom_median)
{
    cmd_arguments args{};
    args.sequence_files = {"./example/mini_example.fasta"};
    args.genome_file = "./example/mini_genom.fasta";
    args.expression_levels = {1};
    args.bits = {1000};
    args.path_out = "./example/";
    args.compressed = true;
    args.k = 4;
    args.window_size = 4;
    args.seed = 0;

    std::vector<uint32_t> expected{4};

    std::vector<uint32_t> medians = ibf(args);

    EXPECT_EQ(expected, medians);
}

TEST(ibf, genom_mean)
{
    cmd_arguments args{};
    args.sequence_files = {"./example/mini_example.fasta"};
    args.genome_file = "./example/mini_genom.fasta";
    args.expression_levels = {1};
    args.bits = {1000};
    args.path_out = "./example/";
    args.aggregate_by = "mean";
    args.compressed = true;
    args.k = 4;
    args.window_size = 4;
    args.seed = 0;

    std::vector<uint32_t> expected{4};

    std::vector<uint32_t> means = ibf(args);

    EXPECT_EQ(expected, means);
}

TEST(search, small_example)
{
    cmd_arguments args{};
    std::vector<uint32_t> expected{1};
    args.sequence_files = {"./example/mini_example.fasta"};
    args.expression_levels = {1};
    args.bits = {1000};
    args.path_out = "./example/";
    args.compressed = true;
    args.k = 4;
    args.window_size = 4;
    args.seed = 0;

    ibf(args);

    args.gene_file = "./example/mini_gen.fasta";
    args.path_in = args.path_out;
    args.expression = 1;

    seqan3::binning_directory_compressed bd;

    std::vector<uint32_t> results{search(bd, args)};

    EXPECT_EQ(expected, results);
}

TEST(search, small_example_gene_not_found)
{
    cmd_arguments args{};
    std::vector<uint32_t> expected{0};
    args.sequence_files = {"./example/mini_example.fasta"};
    args.expression_levels = {1};
    args.bits = {1000};
    args.path_out = "./example/";
    args.compressed = true;
    args.k = 4;
    args.window_size = 4;
    args.seed = 0;

    ibf(args);

    args.gene_file = "./example/mini_gen2.fasta";
    args.path_in = args.path_out;
    args.expression = 1;

    seqan3::binning_directory_compressed bd;

    std::vector<uint32_t> results{search(bd, args)};

    EXPECT_EQ(expected, results);
}

TEST(search, example)
{
    cmd_arguments args{};
    std::vector<uint32_t> expected{0,1};
    // ./needle-ibf ./example/exp_*.fasta -m 2 -m 2 -o ./example/ -e 0.5 -z 1559922 -c
    args.sequence_files = {"./example/exp_01.fasta", "./example/exp_02.fasta", "./example/exp_11.fasta", "./example/exp_12.fasta"};
    args.samples = {2,2};
    args.expression_levels = {0.5};
    args.bits = {1559922};
    args.path_out = "./example/";
    args.compressed = true;
    ibf(args);

    // ./needle-search ./example/gene.fasta -i ./example/ -e 0.5 -c
    args.gene_file = "./example/gene.fasta";
    args.path_in = args.path_out;
    args.expression = 0.5;

    seqan3::binning_directory_compressed bd;

    std::vector<uint32_t> results{search(bd, args)};

    EXPECT_EQ(expected, results);
}
