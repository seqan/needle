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
        "terminate called after throwing an instance of 'seqan3::too_few_arguments'\n"
        "  what():  Please specify which sub program you want to use (one of [ibf,ibfmin,insert,minimiser,search,stats,"
        "test]). Use -h/--help for more information.\nAborted\n"
    };
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}
