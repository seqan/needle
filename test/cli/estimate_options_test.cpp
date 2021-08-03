#include <string>                // strings

#include "cli_test.hpp"

#include "ibf.h"
#include "shared.h"

std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory

TEST_F(cli_test, no_options)
{
    cli_test_result result = execute_app("needle estimate");
    std::string expected
    {
        "needle-estimate - Estimate expression value of transcript based on IBFs.\n"
        "========================================================================\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, fail_no_argument)
{
    cli_test_result result = execute_app("needle estimate", "-m");
    std::string expected
    {
        "Error. Incorrect command line input for estimate. Not enough positional arguments provided "
        "(Need at least 1). See -h/--help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(cli_test, with_argument)
{
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    ibf_args.expression_levels = {1, 2};
    ibf_args.fpr = {0.05};
    std::vector<std::filesystem::path> sequence_files = {std::string(DATA_INPUT_DIR) + "exp_01.fasta"};
    ibf_args.path_out = tmp_dir/"Test_";
    ibf(sequence_files, ibf_args, minimiser_args);

    cli_test_result result = execute_app("needle estimate -i ", tmp_dir/"Test_", data("mini_gen.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});

    std::filesystem::remove(tmp_dir/"Test_IBF_Data");
    std::filesystem::remove(tmp_dir/"Test_IBF_1");
    std::filesystem::remove(tmp_dir/"Test_IBF_2");

}

TEST_F(cli_test, with_argument_normalization_method)
{
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    ibf_args.expression_levels = {1, 2};
    ibf_args.fpr = {0.05};
    std::vector<std::filesystem::path> sequence_files = {std::string(DATA_INPUT_DIR) + "exp_01.fasta"};
    ibf_args.path_out = tmp_dir/"Test_";
    ibf(sequence_files, ibf_args, minimiser_args);

    cli_test_result result = execute_app("needle estimate -m -i ", tmp_dir/"Test_", data("mini_gen.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
    std::filesystem::remove(tmp_dir/"Test_IBF_Data");
    std::filesystem::remove(tmp_dir/"Test_IBF_1");
    std::filesystem::remove(tmp_dir/"Test_IBF_2");
}
