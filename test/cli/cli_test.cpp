#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <cstdlib>               // system calls
#include <seqan3/std/filesystem> // test directory creation
#include <seqan3/std/ranges>     // range comparisons
#include <string>                // strings
#include <sstream>               // ostringstream

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/input.hpp>

// Provides functions for CLI test implementation.
struct cli_test : public ::testing::Test
{
protected:
    // Result struct for returning captured streams and exit code.
    struct cli_test_result
    {
        std::string out{};
        std::string err{};
        int exit_code{};
    };

    // Invoke the app execution. The command must be given in separate parameters.
    template <typename... CommandItemTypes>
    cli_test_result execute_app(CommandItemTypes &&... command_items)
    {
        cli_test_result result{};

        // assemble the command string
        std::ostringstream command{};
        command << "SEQAN3_NO_VERSION_CHECK=1 " << BINDIR;
        int a[] = {0, ((void)(command << command_items << ' '), 0) ... };
        (void) a;

        // always capture the streams
        testing::internal::CaptureStdout();
        testing::internal::CaptureStderr();

        // run the command
        result.exit_code = std::system(command.str().c_str());
        result.out = testing::internal::GetCapturedStdout();
        result.err = testing::internal::GetCapturedStderr();
        return result;
    }

    // Generate the full path of a test file that is provided in the data directory.
    static
    std::filesystem::path data(std::string const & filename)
    {
        return std::filesystem::path{std::string{DATADIR}}.concat(filename);
    }

    // Holds the original work directory.
    std::filesystem::path original_workdir{};

    // Creates a work directory for the test.
    void SetUp() override
    {
        ::testing::TestInfo const * const info = ::testing::UnitTest::GetInstance()->current_test_info();
        std::filesystem::path const test_dir{std::string{OUTPUTDIR} +
                                             std::string{info->test_case_name()} +
                                             std::string{"."} +
                                             std::string{info->name()}};
        try
        {
            std::filesystem::remove_all(test_dir);              // delete the directory if it exists
            std::filesystem::create_directories(test_dir);      // create new empty directory
            original_workdir = std::filesystem::current_path(); // store original wd
            std::filesystem::current_path(test_dir);            // change wd
        }
        catch (std::exception const & exc)
        {
            FAIL() << "Failed to set up the test directory " << test_dir << ":\n" << exc.what();
        }
    }

    // Switches back to the initial working directory.
    void TearDown() override
    {
        try
        {
            std::filesystem::current_path(original_workdir);   // restore original wd
        }
        catch (std::exception const & exc)
        {
            FAIL() << "Failed to set the work directory to " << original_workdir << ":\n" << exc.what();
        }
    }
};

TEST_F(cli_test, no_options)
{
    cli_test_result result = execute_app("fastq_to_fasta");
    std::string expected
    {
        "Fastq-to-Fasta-Converter\n"
        "========================\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, fail_no_argument)
{
    cli_test_result result = execute_app("fastq_to_fasta", "-v");
    std::string expected
    {
        "Parsing error. Not enough positional arguments provided (Need at least 1). "
        "See -h/--help for more information.\n"
    };
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(cli_test, with_argument)
{
    cli_test_result result = execute_app("fastq_to_fasta", data("in.fastq"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, "> seq1\nACGTTTGATTCGCG\n> seq2\nTCGGGGGATTCGCG\n");
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, with_argument_verbose)
{
    cli_test_result result = execute_app("fastq_to_fasta", data("in.fastq"), "-v");
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, "> seq1\nACGTTTGATTCGCG\n> seq2\nTCGGGGGATTCGCG\n");
    EXPECT_EQ(result.err, "Conversion was a success. Congrats!\n");
}

TEST_F(cli_test, with_out_file)
{
    cli_test_result result = execute_app("fastq_to_fasta", data("in.fastq"), "-o", "out.fasta");
    seqan3::sequence_file_input fin{"out.fasta", seqan3::fields<seqan3::field::seq, seqan3::field::id>{}};

    // create records to compare
    using record_type = typename decltype(fin)::record_type;
    using seqan3::operator""_dna5;
    std::vector<record_type> records{};
    records.emplace_back("ACGTTTGATTCGCG"_dna5, std::string{"seq1"});
    records.emplace_back("TCGGGGGATTCGCG"_dna5, std::string{"seq2"});

    EXPECT_TRUE(std::ranges::equal(fin, records));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
}
