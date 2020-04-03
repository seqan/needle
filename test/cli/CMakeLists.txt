cmake_minimum_required (VERSION 3.8)

# The cli_test build target.
add_executable (cli_test cli_test.cpp)
target_link_libraries (cli_test seqan3::seqan3 gmock gmock_main gtest pthread)
add_dependencies(cli_test "${PROJECT_NAME}")
add_test (NAME command_line_interface
          COMMAND cli_test
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
