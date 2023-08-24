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

TEST_F(delete_options_test, with_argument)
{
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    ibf_args.expression_thresholds = {1, 2};
    std::vector<double> fpr = {0.05};
    std::vector<std::filesystem::path> sequence_files = {data("exp_01.fasta"),data("exp_01.fasta")};
    ibf_args.path_out = "Test_";
    std::vector<uint8_t> cutoffs{};
    ibf(sequence_files, ibf_args, minimiser_args, fpr, cutoffs);

    cli_test_result result = execute_app("needle delete -i ", "Test_ ", "0");
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}
