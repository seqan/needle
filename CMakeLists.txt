cmake_minimum_required (VERSION 3.8)

## CUSTOMISE

# Define the application name and version.
project (fastq_to_fasta VERSION 1.0.0)
set (CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Specify the directories where to store the built archives, libraries and executables
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

## BUILD
add_subdirectory(src)
add_executable ("${PROJECT_NAME}" src/main.cpp)
target_link_libraries ("${PROJECT_NAME}" "${PROJECT_NAME}_lib")

## TEST
add_definitions(-DDATADIR=\"${CMAKE_CURRENT_BINARY_DIR}/data/\")
add_definitions(-DBINDIR=\"${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/\")

# Dependency: Googletest
add_subdirectory(lib/gtest)

enable_testing ()
include(test/data/datasources.cmake)
add_subdirectory(test)
