#include <string>                // strings

#include "cli_test.hpp"

struct count_options_test : public cli_test {};

TEST_F(count_options_test, no_options)
{
    cli_test_result result = execute_app("needle count");
    std::string expected
    {
        "needle-count - Get expression value depending on minimizers. This function is only "
        "used to test the validity of Needle's estimation approach. It "
        "estimates the expression value for all sequences in the genome file "
        "based on the exact minimiser occurrences of the given sequence files.\n"
        "======================================================================"
        "======================================================================"
        "======================================================================"
        "======================================================================"
        "==\n    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(count_options_test, fail_no_argument)
{
    cli_test_result result = execute_app("needle count", "--seed 0");
    std::string expected
    {
        "Error. Incorrect command line input for count. Not enough positional arguments provided "
        "(Need at least 1). See -h/--help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(count_options_test, with_arguments)
{
    cli_test_result result = execute_app("needle count -k 4 -w 4 --seed 0 -g", data("mini_gen.fasta"), data("mini_example.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(count_options_test, multithreads)
{
    cli_test_result result = execute_app("needle count -k 4 -w 8 -t 2 -g",
                                        data("mini_gen.fasta"),
                                        data("mini_example.fasta"),
                                        data("mini_example.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(count_options_test, paired)
{
    cli_test_result result = execute_app("needle count -k 4 -w 8 -p -g",
                                        data("mini_gen.fasta"),
                                        data("mini_example.fasta"),
                                        data("mini_example.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}
