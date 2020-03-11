cmake_minimum_required (VERSION 3.4)

# The api_test build target.
add_executable (test-needle test.cpp)
target_link_libraries(test-needle gtest_main)
target_include_directories (test-needle PUBLIC ../../include)
target_link_libraries (test-needle seqan3::seqan3)
target_link_libraries (test-needle robin_hood)
add_test (test-needle test-needle)
