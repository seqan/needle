#include <string>                // strings

#include "cli_test.hpp"

TEST_F(cli_test, no_options)
{
    cli_test_result result = execute_app("needle insert");
    std::string expected
    {
        "needle-insert - Inserts an IBF to another IBF.\n"
        "==============================================\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, fail_no_argument)
{
    cli_test_result result = execute_app("needle insert", "-c");
    std::string expected
    {
        "Error. Incorrect command line input for insert. Not enough positional arguments provided "
        "(Need at least 1). See -h/--help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}
