#include <string>                // strings

#include "cli_test.hpp"

#include "ibf.h"
#include "shared.h"

struct delete_options_test : public cli_test {};

TEST_F(delete_options_test, delete_no_options)
{
    cli_test_result result = execute_app("needle delete");
    std::string expected
    {
        "needle-delete - Delete experiments specified by their position from the Needle index.\n"
        "=====================================================================================\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}
