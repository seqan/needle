#include <string>                // strings

#include "cli_test.hpp"

TEST_F(cli_test, no_options)
{
    cli_test_result result = execute_app("needle stats");
    std::string expected
    {
        "needle-stats - Get statistics from header files produced by needle minimiser.\n"
        "=============================================================================\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, with_argument)
{
    cli_test_result result = execute_app("needle stats", data("mini_example.header"), data( "mini_example2.header"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.err, std::string{});
}
