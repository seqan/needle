cmake_minimum_required (VERSION 3.4)

# The api_test build target.
add_executable (fastq_conversion_test fastq_conversion_test.cpp)
target_link_libraries (fastq_conversion_test "${PROJECT_NAME}_object" gmock gmock_main gtest pthread)
add_test (fastq_conversion_test fastq_conversion_test)
