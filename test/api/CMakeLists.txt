cmake_minimum_required (VERSION 3.8)

# The api_test build target.
add_executable (api_test fastq_conversion_test.cpp)
target_link_libraries (api_test "${PROJECT_NAME}_lib" gmock gmock_main gtest pthread)
add_test (NAME unit COMMAND api_test)
