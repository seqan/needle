#include <string>                // strings

#include "../app_test.hpp"

struct ibf_options_test : public app_test {};

TEST_F(ibf_options_test, ibf_no_options)
{
    app_test_result result = execute_app("ibf");
    std::string expected
    {
        "needle-ibf - Constructs the Needle index.\n"
        "=========================================\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(ibf_options_test, ibf_fail_no_argument)
{
    app_test_result result = execute_app("ibf", "-c");
    std::string expected
    {
        "Error. Incorrect command line input for ibf. Not enough positional arguments provided "
        "(Need at least 1). See -h/--help for more information.\n"
    };
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, std::string{""});
    EXPECT_EQ(result.err, expected);
}

TEST_F(ibf_options_test, ibf_fail_contradiction)
{
    app_test_result result = execute_app("ibf -f 0.05 -e 1 -e 2 -l 1", data("exp_01.fasta"));
    std::string expected
    {
        "Error. Please set the expression levels OR give the number of expression levels.\n"
    };
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, std::string{""});
    EXPECT_EQ(result.err, expected);
}

TEST_F(ibf_options_test, ibf_fail_contradiction2)
{
    app_test_result result = execute_app("ibf -f 0.05 -e 1 -e 2 --levels-by-genome ", data("exp_01.fasta"),
                                         data("exp_01.fasta"));
    std::string expected
    {
        "Error. The determination of expression levels can not be used with individual levels already given. Please set "
        "the expression levels without the option --level-by-genome OR use the number of expression levels with that option."
        "\n"
    };
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, std::string{""});
    EXPECT_EQ(result.err, expected);
}

TEST_F(ibf_options_test, ibf_fail_no_fpr)
{
    app_test_result result = execute_app("ibf -l 2", data("exp_01.fasta"));
    std::string expected
    {
        "Error. Please give a false positive rate for the IBFs.\n"
    };
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, std::string{""});
    EXPECT_EQ(result.err, expected);
}

TEST_F(ibf_options_test, ibf_fail_incorrect_number_of_fprs)
{
    app_test_result result = execute_app("ibf -f 0.05 -f 0.01 -f 0.03 -l 2", data("exp_01.fasta"));
    std::string expected
    {
        "Error. Length of false positive rates for IBFs is not equal to length of expression thresholds.\n"
    };
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, std::string{""});
    EXPECT_EQ(result.err, expected);
}

TEST_F(ibf_options_test, ibf_with_argument)
{
    app_test_result result = execute_app("ibf -f 0.05 -l 1", data("exp_01.fasta"));
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(ibf_options_test, ibf_with_argument_with_options)
{
    app_test_result result = execute_app("ibf -f 0.05 -k 4 -w 8 -l 1", data("exp_01.fasta"));
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(ibf_options_test, ibf_with_argument_with_userdefined_shape)
{
    app_test_result result = execute_app("ibf -f 0.05 -k 4 -w 8 --shape 11 -l 1", data("exp_01.fasta"));
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(ibf_options_test, ibfmin_no_options)
{
    app_test_result result = execute_app("ibfmin");
    std::string expected
    {
        "needle-ibfmin - Constructs the Needle index from the minimiser files created by needle minimiser.\n"
        "=================================================================================================\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(ibf_options_test, ibfmin_fail_no_argument)
{
    app_test_result result = execute_app("ibfmin -c");
    std::string expected
    {
        "Error. Incorrect command line input for ibfmin. Not enough positional arguments provided "
        "(Need at least 1). See -h/--help for more information.\n"
    };
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(ibf_options_test, ibfmin_fail_contradiction)
{
    app_test_result result = execute_app("ibfmin -f 0.05 -e 1 -e 2 -l 1", data("mini_example.minimiser"));
    std::string expected
    {
        "Error. Please set the expression levels OR give the number of expression levels.\n"
    };
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, std::string{""});
    EXPECT_EQ(result.err, expected);
}

TEST_F(ibf_options_test, ibfmin_fail_contradiction2)
{
    app_test_result result = execute_app("ibfmin -f 0.05 -e 1 -e 2 --levels-by-genome ", data("exp_01.fasta"),
                                         data("mini_example.minimiser"));
    std::string expected
    {
        "Error. The determination of expression levels can not be used with individual levels already given. Please set "
        "the expression levels without the option --level-by-genome OR use the number of expression levels with that option."
        "\n"
    };
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, std::string{""});
    EXPECT_EQ(result.err, expected);
}

TEST_F(ibf_options_test, ibfmin_fail_no_fpr)
{
    app_test_result result = execute_app("ibfmin -l 2", data("mini_example.minimiser"));
    std::string expected
    {
        "Error. Please give a false positive rate for the IBFs.\n"
    };
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, std::string{""});
    EXPECT_EQ(result.err, expected);
}

TEST_F(ibf_options_test, ibfmin_fail_incorrect_number_of_fprs)
{
    app_test_result result = execute_app("ibfmin -f 0.05 -f 0.01 -f 0.03 -l 2", data("mini_example.minimiser"));
    std::string expected
    {
        "Error. Length of false positive rates for IBFs is not equal to length of expression thresholds.\n"
    };
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, std::string{""});
    EXPECT_EQ(result.err, expected);
}

TEST_F(ibf_options_test, ibfmin_with_argument)
{
    app_test_result result = execute_app("ibfmin -f 0.05 -l 1", data("mini_example.minimiser"));
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(ibf_options_test, compressed)
{
    app_test_result result = execute_app("ibfmin -f 0.05 -l 1 -c ", data("mini_example.minimiser"));
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(ibf_options_test, more_hash_functions)
{
    app_test_result result = execute_app("ibfmin -f 0.05 -l 1 -n 4 ", data("mini_example.minimiser"));
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(ibf_options_test, expression_thresholds)
{
    app_test_result result = execute_app("ibfmin -f 0.05 -e 2 -e 4", data("mini_example.minimiser"));
    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, std::string{});
}
