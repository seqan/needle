#include <string> // strings

#include "../app_test.hpp"

struct needle_options_test : public app_test
{};

TEST_F(needle_options_test, no_options)
{
    app_test_result result = execute_app();
    std::string expected{"needle\n"
                         "======\n"
                         "    Try -h or --help for more information.\n"};
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(needle_options_test, fail_no_argument)
{
    app_test_result result = execute_app("-v");
    std::string expected{
        "Error. Incorrect command. See needle help for more information.You either forgot or misspelled the subcommand!"
        " Please specify which sub-program you want to use: one of "
        "[count,delete,estimate,genome,ibf,ibfmin,insert,insertmin,minimiser]. "
        "Use -h/--help for more information.\n"};
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_TRUE(result.err == expected);
}
