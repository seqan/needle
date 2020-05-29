#include <string>                // strings

#include "cli_test.hpp"

TEST_F(cli_test, no_options)
{
    cli_test_result result = execute_app("needle test");
    std::string expected
    {
        "needle-test - Constructs multiple IBFs for given experiments using different normalization methods.\n"
        "===================================================================================================\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, fail_no_argument)
{
    cli_test_result result = execute_app("needle test", "-c");
    std::string expected
    {
        "Error. Incorrect command line input for test. Not enough positional arguments provided "
        "(Need at least 1). See -h/--help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(cli_test, with_argument)
{
    cli_test_result result = execute_app("needle test -k 4 -w 8 -g",  data("mini_genom.fasta"), data("mini_example.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}
