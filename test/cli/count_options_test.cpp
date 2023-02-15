#include <string>                // strings

#include "cli_test.hpp"

struct count_options_test : public cli_test {};

TEST_F(count_options_test, no_options)
{
    cli_test_result result = execute_app("needle count");
    std::string expected
    {
        "needle-count - Get expression value depending on minimizers. "
        "This function is an alternative to pseudoaligners like kallisto. It "
        "estimates the expression value for all sequences in the genome file "
        "based on the exact minimiser occurrences of the given sequence files. "
        "Please run genome beforehand to create the genome file.\n"
        "======================================================================"
        "======================================================================"
        "======================================================================"
        "======================================================================"
        "==========================================\n    "
        "Try -h or --help for more information.\n"
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
    cli_test_result result = execute_app("needle count -k 4 -w 4 --seed 0 --genome", data("mini_gen.genome")," --include", data("mini_gen.fasta"), data("mini_example.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(count_options_test, multithreads)
{
    cli_test_result result = execute_app("needle count -k 4 -w 8 -t 2 --genome",
                                        data("mini_gen.genome"),
                                        " --include",
                                        data("mini_gen.fasta"),
                                        data("mini_example.fasta"),
                                        data("mini_example.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(count_options_test, paired)
{
    cli_test_result result = execute_app("needle count -k 4 -w 8 -p --genome",
                                        data("mini_gen.genome"),
                                        " --include",
                                        data("mini_gen.fasta"),
                                        data("mini_example.fasta"),
                                        data("mini_example.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}

struct genome_options_test : public cli_test {};

TEST_F(genome_options_test, no_options)
{
    cli_test_result result = execute_app("needle genome");
    std::string expected
    {
        "needle-genome - Creates the genome file necessary as an input to count.\n"
        "=======================================================================\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(genome_options_test, fail_no_argument)
{
    cli_test_result result = execute_app("needle genome", "--seed 0");
    std::string expected
    {
        "Error. Incorrect command line input for count. Not enough positional arguments provided "
        "(Need at least 1). See -h/--help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(genome_options_test, with_arguments)
{
    cli_test_result result = execute_app("needle genome -k 4 -w 4 --seed 0 ", data("mini_gen.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}
