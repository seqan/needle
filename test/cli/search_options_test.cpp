#include <string>                // strings

#include "cli_test.hpp"

TEST_F(cli_test, no_options)
{
    cli_test_result result = execute_app("needle search");
    std::string expected
    {
        "needle-search - Search through an IBF.\n"
        "======================================\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, fail_no_argument)
{
    cli_test_result result = execute_app("needle search", "-c");
    std::string expected
    {
        "Error. Incorrect command line input for search. Not enough positional arguments provided "
        "(Need at least 1). See -h/--help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(cli_test, with_argument)
{
    cli_test_result result = execute_app("needle search -k 4 -w 4 -s 0 -e 1 -c -i ", data(""), data("mini_gen.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, "Results:\n1 \n");
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, with_argument_gene_not_found)
{
    cli_test_result result = execute_app("needle search -k 4 -w 4 -s 0 -e 1 -c -i ", data(""), data("mini_gen2.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, "Results:\n0 \n");
    EXPECT_EQ(result.err, std::string{});
}

// https://github.com/MitraDarja/needle/issues/33
TEST_F(cli_test, issue33)
{
    cli_test_result result = execute_app("needle search -k 4 -w 4 -s 0 -t 0.7 -e 1 -c -i ", data(""),
                                         data("mini_gen.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, "Results:\n1 \n");
    EXPECT_EQ(result.err, std::string{});
}
