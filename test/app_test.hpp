// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#pragma once

#include <gtest/gtest.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

// Checks for CLI test result for success, and prints the command line call if the test fails.
#ifndef EXPECT_SUCCESS
#    define EXPECT_SUCCESS(arg)                                                                                        \
        EXPECT_EQ(arg.exit_code, 0) << "Command: " << arg.command << "\n Working directory: " << arg.current_workdir
#endif

// Checks the exit code of a CLI test result for failure, and prints the command line call if the test fails.
#ifndef EXPECT_FAILURE
#    define EXPECT_FAILURE(arg)                                                                                        \
        EXPECT_NE(arg.exit_code, 0) << "Command: " << arg.command << "\n Working directory: " << arg.current_workdir
#endif

// Provides functions for CLI test implementation.
struct app_test : public ::testing::Test
{
private:
    // Holds the original work directory where Gtest has been started.
    std::filesystem::path original_workdir{};

    // Holds the current work directory.
    std::filesystem::path current_workdir{};

protected:
    // Result struct for captured streams and exit code.
    struct app_test_result
    {
        std::string out{};
        std::string err{};
        std::string command{};
        std::string current_workdir{};
        int exit_code{};
    };

    // Invoke the app execution. The command line call should be given as separate parameters.
    template <typename... CommandItemTypes>
    app_test_result execute_app(CommandItemTypes &&... command_items)
    {
        app_test_result result{.current_workdir = current_workdir};

        // Assemble the command string and disable version check.
        result.command = [&command_items...]()
        {
            std::ostringstream command{};
            command << "SHARG_NO_VERSION_CHECK=1 " << BINDIR << APPNAME;
            (void)((command << ' ' << command_items), ...); // (void) silences "state has no effect" warning if no items
            return command.str();
        }();

        // Always capture the output streams.
        testing::internal::CaptureStdout();
        testing::internal::CaptureStderr();

        // Run the command and return results.
        result.exit_code = std::system(result.command.c_str());
        result.out = testing::internal::GetCapturedStdout();
        result.err = testing::internal::GetCapturedStderr();
        return result;
    }

    // Generate the full path of a test input file that is provided in the data directory.
    static std::filesystem::path data(std::string const & filename)
    {
        return std::filesystem::path{std::string{DATADIR}}.concat(filename);
    }

    // Read the contents of a file into a string.
    static std::string const string_from_file(std::filesystem::path const & path,
                                              std::ios_base::openmode const mode = std::ios_base::in)
    {
        std::ifstream file_stream(path, mode);
        if (!file_stream.is_open())
            throw std::logic_error{"Cannot open " + path.string()};
        std::stringstream file_buffer;
        file_buffer << file_stream.rdbuf();
        return {file_buffer.str()};
    }

    // Create an individual work directory for the current test.
    void SetUp() override
    {
        // Assemble the directory name.
        ::testing::TestInfo const * const info = ::testing::UnitTest::GetInstance()->current_test_info();
        current_workdir = std::filesystem::path{std::string{OUTPUTDIR} + std::string{info->test_case_name()}
                                                + std::string{"."} + std::string{info->name()}};
        try
        {
            std::filesystem::remove_all(current_workdir);         // delete the directory if it exists
            std::filesystem::create_directories(current_workdir); // create the new empty directory
            original_workdir = std::filesystem::current_path();   // store original work dir path
            std::filesystem::current_path(current_workdir);       // change the work dir
        }
        catch (std::exception const & exc)
        {
            FAIL() << "Failed to set up the test directory " << current_workdir << ":\n" << exc.what();
        }
    }

    // Switch back to the initial work directory.
    void TearDown() override
    {
        try
        {
            std::filesystem::current_path(original_workdir); // restore the original work dir
        }
        catch (std::exception const & exc)
        {
            FAIL() << "Failed to set the work directory to " << original_workdir << ":\n" << exc.what();
        }
    }
};
