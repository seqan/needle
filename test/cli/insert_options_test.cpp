#include <string>                // strings

#include "../app_test.hpp"

#include "ibf.h"
#include "shared.h"

struct insert_options_test : public app_test {};

TEST_F(insert_options_test, insert_no_options)
{
    app_test_result result = execute_app("insert");
    std::string expected
    {
        "needle-insert - Inserts into a given uncompressed Needle index.\n"
        "===============================================================\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(insert_options_test, insert_fail_no_argument)
{
    app_test_result result = execute_app("insert", "-c");
    std::string expected
    {
        "Error. Incorrect command line input for insert. Not enough positional arguments provided "
        "(Need at least 1). See -h/--help for more information.\n"
    };
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, std::string{""});
    EXPECT_EQ(result.err, expected);
}

TEST_F(insert_options_test, with_argument)
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path();
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    ibf_args.expression_thresholds = {1, 2};
    std::vector<double> fpr = {0.05};
    std::vector<std::filesystem::path> sequence_files = {data("exp_01.fasta")};
    ibf_args.path_out = tmp_dir/"Test_";
    std::vector<uint8_t> cutoffs{1};
    ibf(sequence_files, ibf_args, minimiser_args, fpr, cutoffs);

    app_test_result result = execute_app("insert -i ", tmp_dir/"Test_", data("exp_01.fasta"));
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}
