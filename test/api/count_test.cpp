// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <iostream>

#include <seqan3/test/expect_range_eq.hpp>

#include "../app_test.hpp"
#include "ibf.hpp"
#include "shared.hpp"

// To prevent issues when running multiple CLI tests in parallel, give each CLI test unique names:
struct count_test : public app_test
{
    void initialization_args(estimate_ibf_arguments & args)
    {
        args.compressed = true;
        args.k = 4;
        args.shape = seqan3::ungapped{args.k};
        args.w_size = seqan3::window_size{4};
        args.s = seqan3::seed{0};
        args.path_out = "Count_Test_";
    }
};

TEST_F(count_test, small_example)
{
    estimate_ibf_arguments args{};
    initialization_args(args);

    count(args, {data("mini_example.fasta")}, data("mini_gen.fasta"), data("mini_gen.genome"), false);

    std::ifstream output_file("mini_example.count.out");
    std::string line;
    std::string expected{"gen1\t3"};
    if (output_file.is_open())
    {
        while (std::getline(output_file, line))
        {
            EXPECT_EQ(expected, line);
        }
        output_file.close();
    }
}

TEST_F(count_test, small_example_paired)
{
    estimate_ibf_arguments args{};
    initialization_args(args);

    count(args,
          {data("mini_example.fasta"), data("mini_example.fasta")},
          data("mini_gen.fasta"),
          data("mini_gen.genome"),
          true);

    std::ifstream output_file("mini_example.count.out");
    std::string line;
    std::string expected{"gen1\t6"};
    if (output_file.is_open())
    {
        while (std::getline(output_file, line))
        {
            EXPECT_EQ(expected, line);
        }
        output_file.close();
    }
}

TEST_F(count_test, small_example_exclude)
{
    estimate_ibf_arguments args{};
    initialization_args(args);

    count(args, {data("mini_example.fasta")}, data("mini_gen.fasta"), data("mini_gen2.genome"), false);

    std::ifstream output_file("mini_example.count.out");
    std::string line;
    std::string expected{"gen1\t3"};
    if (output_file.is_open())
    {
        while (std::getline(output_file, line))
        {
            EXPECT_EQ(expected, line);
        }
        output_file.close();
    }
}

TEST_F(count_test, genome_small_example)
{
    estimate_ibf_arguments args{};
    initialization_args(args);

    count_genome(args, data("mini_gen.fasta"), data("mini_gen2.fasta"));

    std::ifstream output_file;
    uint64_t expected{192};
    output_file.open("mini_gen.genome", std::ios::binary);
    uint64_t minimiser;
    while (output_file.read((char *)&minimiser, sizeof(minimiser)))
    {
        EXPECT_EQ(expected, minimiser);
    }
    output_file.close();
}
