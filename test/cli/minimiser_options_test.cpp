#include <string>                // strings

#include "cli_test.hpp"

TEST_F(cli_test, no_options)
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

TEST_F(cli_test, fail_no_argument)
{
    cli_test_result result = execute_app("needle minimiser", "-s 0");
    std::string expected
    {
        "Error. Incorrect command line input for minimiser. Not enough positional arguments provided "
        "(Need at least 1). See -h/--help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(cli_test, with_arguments)
{
    cli_test_result result = execute_app("needle minimiser -k 4 -w 4 -s 0", data("mini_example.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, cutoff)
{
    cli_test_result result = execute_app("needle minimiser -k 4 -w 8 -s 0 -u 2", data("mini_example.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, multiple_sample)
{
    cli_test_result result = execute_app("needle minimiser -k 4 -w 8 -s 0 -m 2 ", data("mini_example.fasta"),
                                                                                  data("mini_example.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, multithreads)
{
    cli_test_result result = execute_app("needle minimiser -k 4 -w 8 -s 0 -t 2", data("mini_example.fasta"),
                                                                                 data("mini_example.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, paired)
{
    cli_test_result result = execute_app("needle minimiser -k 4 -w 8 -s 0 -q ", data("mini_example.fasta"),
                                                                                data("mini_example.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}
