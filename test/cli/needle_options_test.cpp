#include <string>                // strings

#include "cli_test.hpp"

TEST_F(cli_test, no_options)
{
    cli_test_result result = execute_app("needle");
    std::string expected
    {
        "needle\n"
        "======\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, fail_no_argument)
{
    cli_test_result result = execute_app("needle", "-v");
    std::string expected
    {
        "Error. Incorrect command. See needle help for more information.You either forgot or misspelled the subcommand!"
        " Please specify which sub-program you want to use: one of [count,estimate,ibf,ibfmin,minimiser]. "
        "Use -h/--help for more information.\n"
    };
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_TRUE(result.err == expected);
}
