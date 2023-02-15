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

void initialization_args(estimate_ibf_arguments & args)
{
    args.compressed = true;
    args.k = 4;
    args.shape = seqan3::ungapped{args.k};
    args.w_size = seqan3::window_size{4};
    args.s = seqan3::seed{0};
    args.path_out = tmp_dir/"Count_Test_";
}

TEST(count, small_example)
{
    estimate_ibf_arguments args{};
    initialization_args(args);

    count(args, {std::string(DATA_INPUT_DIR) + "mini_example.fasta"}, std::string(DATA_INPUT_DIR) + "mini_gen.fasta",
          std::string(DATA_INPUT_DIR) + "mini_gen.genome", false);

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
    std::filesystem::remove(tmp_dir/"Count_Test_mini_example.count.out");
}

TEST(count, small_example_paired)
{
    estimate_ibf_arguments args{};
    initialization_args(args);

    count(args, {std::string(DATA_INPUT_DIR) + "mini_example.fasta", std::string(DATA_INPUT_DIR) + "mini_example.fasta"},
          std::string(DATA_INPUT_DIR) + "mini_gen.fasta", std::string(DATA_INPUT_DIR) + "mini_gen.genome", true);

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
    std::filesystem::remove(tmp_dir/"Count_Test_mini_example.count.out");
}

TEST(count, small_example_exclude)
{
    estimate_ibf_arguments args{};
    initialization_args(args);

    count(args, {std::string(DATA_INPUT_DIR) + "mini_example.fasta"}, std::string(DATA_INPUT_DIR) + "mini_gen.fasta",
                 std::string(DATA_INPUT_DIR) + "mini_gen2.genome", false);

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
    std::filesystem::remove(tmp_dir/"Count_Test_mini_example.count.out");
}

TEST(genome, small_example)
{
    estimate_ibf_arguments args{};
    initialization_args(args);

    count_genome(args, std::string(DATA_INPUT_DIR) + "mini_gen.fasta", std::string(DATA_INPUT_DIR) + "mini_gen2.fasta");

    std::ifstream output_file;
    uint64_t expected{192};
    output_file.open(tmp_dir/"mini_gen.genome", std::ios::binary);
    uint64_t minimiser;
    while(output_file.read((char*)&minimiser, sizeof(minimiser)))
    {
        EXPECT_EQ(expected, minimiser);
    }
    output_file.close();
}
