#include <string>                // strings

#include "cli_test.hpp"

struct minimiser_options_test : public cli_test {};

TEST_F(minimiser_options_test, no_options)
{
    cli_test_result result = execute_app("needle minimiser");
    std::string expected
    {
        "needle-minimiser - Calculates minimiser for given experiments.\n"
        "==============================================================\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(minimiser_options_test, fail_no_argument)
{
    cli_test_result result = execute_app("needle minimiser", "--seed 0");
    std::string expected
    {
        "Error. Incorrect command line input for minimiser. Not enough positional arguments provided "
        "(Need at least 1). See -h/--help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(minimiser_options_test, with_arguments)
{
    cli_test_result result = execute_app("needle minimiser -k 4 -w 4 --seed 0", data("mini_example.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(minimiser_options_test, cutoff)
{
    cli_test_result result = execute_app("needle minimiser -k 4 -w 8 --cutoff 2", data("mini_example.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(minimiser_options_test, multiple_sample)
{
    cli_test_result result = execute_app("needle minimiser -k 4 -w 8 --samples 2 ", data("mini_example.fasta"),
                                                                                   data("mini_example.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(minimiser_options_test, multithreads)
{
    cli_test_result result = execute_app("needle minimiser -k 4 -w 8 -t 2", data("mini_example.fasta"),
                                                                            data("mini_example.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(minimiser_options_test, paired)
{
    cli_test_result result = execute_app("needle minimiser -k 4 -w 8 -p ", data("mini_example.fasta"),
                                                                           data("mini_example.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}
