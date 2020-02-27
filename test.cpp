#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>

#include "ibf.h"
#include "minimizer.h"
#include "search.h"

void initialization_args(arguments & args)
{
    args.compressed = true;
    args.k = 4;
    args.window_size = 4;
    args.seed = 0;
}

void initialization_ibf_args(ibf_arguments & args)
{
    args.expression_levels = {1};
    args.bits = {1000};
    args.path_out = "./example/";
}

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
    arguments args{};
    ibf_arguments ibf_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    ibf_args.sequence_files = {"./example/mini_example.fasta"};

    std::vector<uint32_t> expected{3};

    std::vector<uint32_t> medians = ibf(args, ibf_args);

    EXPECT_EQ(expected, medians);
}

TEST(ibf, mean)
{
    arguments args{};
    ibf_arguments ibf_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    ibf_args.sequence_files = {"./example/mini_example.fasta"};
    ibf_args.aggregate_by = "mean";

    std::vector<uint32_t> expected{3};

    std::vector<uint32_t> means = ibf(args, ibf_args);

    EXPECT_EQ(expected, means);
}

TEST(ibf, random)
{
    arguments args{};
    ibf_arguments ibf_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    ibf_args.sequence_files = {"./example/mini_example.fasta"};
    ibf_args.aggregate_by = "random";
    ibf_args.random = 40;

    std::vector<uint32_t> medians = ibf(args, ibf_args);

    EXPECT_TRUE((3 == medians[0]));
}

TEST(ibf, genom_median)
{
    arguments args{};
    ibf_arguments ibf_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    ibf_args.sequence_files = {"./example/mini_example.fasta"};
    ibf_args.genome_file = "./example/mini_genom.fasta";

    std::vector<uint32_t> expected{4};

    std::vector<uint32_t> medians = ibf(args, ibf_args);

    EXPECT_EQ(expected, medians);
}

// Test a genome which has so different minimizers that some minimizers in seq can not be found
TEST(ibf, genom_median_no_match)
{
    arguments args{};
    ibf_arguments ibf_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    ibf_args.sequence_files = {"./example/mini_example2.fasta"};
    ibf_args.genome_file = "./example/mini_genom.fasta";

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
    ibf_args.sequence_files = {"./example/mini_example.fasta"};
    ibf_args.genome_file = "./example/mini_genom.fasta";
    ibf_args.aggregate_by = "mean";

    std::vector<uint32_t> expected{4};

    std::vector<uint32_t> means = ibf(args, ibf_args);

    EXPECT_EQ(expected, means);
}

TEST(search, small_example)
{
    arguments args{};
    ibf_arguments ibf_args{};
    search_arguments search_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    std::vector<uint32_t> expected{1};
    ibf_args.sequence_files = {"./example/mini_example.fasta"};

    ibf(args, ibf_args);

    search_args.gene_file = "./example/mini_gen.fasta";
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
    ibf_args.sequence_files = {"./example/mini_example.fasta"};

    ibf(args, ibf_args);

    search_args.gene_file = "./example/mini_gen2.fasta";
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
    ibf_args.sequence_files = {"./example/mini_example.fasta"};
    ibf_args.expression_levels = {0};
    ibf_args.cutoffs = {2};

    ibf(args, ibf_args);

    search_args.gene_file = "./example/mini_gen3.fasta";
    search_args.path_in = ibf_args.path_out;
    search_args.expression = 1;

    std::vector<uint32_t> results{search(args, search_args)};

    EXPECT_EQ(expected, results);
}

TEST(search, example)
{
    arguments args{};
    ibf_arguments ibf_args{};
    search_arguments search_args{};
    std::vector<uint32_t> expected{0,1};
    // ./needle-ibf ./example/exp_*.fasta -m 2 -m 2 -o ./example/ -e 0.5 -z 1559922 -c
    ibf_args.sequence_files = {"./example/exp_01.fasta", "./example/exp_02.fasta", "./example/exp_11.fasta", "./example/exp_12.fasta"};
    ibf_args.samples = {2,2};
    ibf_args.expression_levels = {0.5};
    ibf_args.bits = {1559922};
    ibf_args.path_out = "./example/";
    args.compressed = true;
    ibf(args, ibf_args);

    // ./needle-search ./example/gene.fasta -i ./example/ -e 0.5 -c
    search_args.gene_file = "./example/gene.fasta";
    search_args.path_in = ibf_args.path_out;
    search_args.expression = 0.5;

    std::vector<uint32_t> results{search(args, search_args)};

    EXPECT_EQ(expected, results);
}
