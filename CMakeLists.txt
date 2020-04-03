cmake_minimum_required (VERSION 3.8)

## CUSTOMISE

# Define the application name and version.
project (fastq_to_fasta VERSION 1.0.0)
set (CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


## BUILD
add_subdirectory(src)
add_executable ("${PROJECT_NAME}" src/main.cpp)
target_link_libraries ("${PROJECT_NAME}" "${PROJECT_NAME}_lib")


## TEST

# copy test data and define DATADIR
file(COPY test/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/test)
add_definitions(-DDATADIR=\"${CMAKE_CURRENT_BINARY_DIR}/test/data/\")
add_definitions(-DBINDIR=\"${CMAKE_CURRENT_BINARY_DIR}/\")

# Dependency: Googletest
add_subdirectory(lib/gtest)

enable_testing ()
add_subdirectory(test/api)
add_subdirectory(test/cli)
