#include <string>                // strings

#include "cli_test.hpp"

struct ibf_options_test : public cli_test {};

TEST_F(ibf_options_test, ibf_no_options)
{
    cli_test_result result = execute_app("needle ibf");
    std::string expected
    {
        "needle-ibf - Constructs the Needle index.\n"
        "=========================================\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(ibf_options_test, ibf_fail_no_argument)
{
    cli_test_result result = execute_app("needle ibf", "-c");
    std::string expected
    {
        "Error. Incorrect command line input for ibf. Not enough positional arguments provided "
        "(Need at least 1). See -h/--help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{""});
    EXPECT_EQ(result.err, expected);
}

TEST_F(ibf_options_test, ibf_fail_contradiction)
{
    cli_test_result result = execute_app("needle ibf -f 0.05 -e 1 -e 2 -l 1", data("exp_01.fasta"));
    std::string expected
    {
        "Error. Please set the expression levels OR give the number of expression levels.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{""});
    EXPECT_EQ(result.err, expected);
}

TEST_F(ibf_options_test, ibf_fail_contradiction2)
{
    cli_test_result result = execute_app("needle ibf -f 0.05 -e 1 -e 2 --levels-by-genome ", data("exp_01.fasta"),
                                         data("exp_01.fasta"));
    std::string expected
    {
        "Error. The determination of expression levels can not be used with individual levels already given. Please set "
        "the expression levels without the option --level-by-genome OR use the number of expression levels with that option."
        "\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{""});
    EXPECT_EQ(result.err, expected);
}

TEST_F(ibf_options_test, ibf_with_argument)
{
    cli_test_result result = execute_app("needle ibf -f 0.05 -l 1", data("exp_01.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(ibf_options_test, ibf_with_argument_with_options)
{
    cli_test_result result = execute_app("needle ibf -f 0.05 -k 4 -w 8 -l 1", data("exp_01.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(ibf_options_test, ibf_with_argument_with_userdefined_shape)
{
    cli_test_result result = execute_app("needle ibf -f 0.05 -k 4 -w 8 --shape 2 -l 1", data("exp_01.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(ibf_options_test, ibfmin_no_options)
{
    cli_test_result result = execute_app("needle ibfmin");
    std::string expected
    {
        "needle-ibfmin - Constructs the Needle index from the minimiser files created by needle minimiser.\n"
        "=================================================================================================\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(ibf_options_test, ibfmin_fail_no_argument)
{
    cli_test_result result = execute_app("needle ibfmin -c");
    std::string expected
    {
        "Error. Incorrect command line input for ibfmin. Not enough positional arguments provided "
        "(Need at least 1). See -h/--help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(ibf_options_test, ibfmin_fail_contradiction)
{
    cli_test_result result = execute_app("needle ibfmin -f 0.05 -e 1 -e 2 -l 1", data("mini_example.minimiser"));
    std::string expected
    {
        "Error. Please set the expression levels OR give the number of expression levels.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{""});
    EXPECT_EQ(result.err, expected);
}

TEST_F(ibf_options_test, ibfmin_fail_contradiction2)
{
    cli_test_result result = execute_app("needle ibfmin -f 0.05 -e 1 -e 2 --levels-by-genome ", data("exp_01.fasta"),
                                         data("mini_example.minimiser"));
    std::string expected
    {
        "Error. The determination of expression levels can not be used with individual levels already given. Please set "
        "the expression levels without the option --level-by-genome OR use the number of expression levels with that option."
        "\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{""});
    EXPECT_EQ(result.err, expected);
}

TEST_F(ibf_options_test, ibfmin_with_argument)
{
    cli_test_result result = execute_app("needle ibfmin -f 0.05 -l 1", data("mini_example.minimiser"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(ibf_options_test, compressed)
{
    cli_test_result result = execute_app("needle ibfmin -f 0.05 -l 1 -c ", data("mini_example.minimiser"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(ibf_options_test, more_hash_functions)
{
    cli_test_result result = execute_app("needle ibfmin -f 0.05 -l 1 -n 4 ", data("mini_example.minimiser"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(ibf_options_test, expression_thresholds)
{
    cli_test_result result = execute_app("needle ibfmin -f 0.05 -e 2 -e 4", data("mini_example.minimiser"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}
