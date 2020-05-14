cmake_minimum_required (VERSION 3.8)

## CUSTOMISE

# Define the application name and version.
project (needle VERSION 1.0.0)
set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

## BUILD

# Specify the directories where to store the built archives, libraries and executables
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
add_subdirectory(src)
add_subdirectory(lib/robin-hood-hashing)
set(DATA_DIR DATADIR=${CMAKE_CURRENT_SOURCE_DIR})
add_compile_definitions(DATA_DIR="${CMAKE_CURRENT_LIST_DIR}/test/data/")
set(CMAKE_BUILD_TYPE Release)
add_executable (needle src/needle.cpp)
target_link_libraries (needle seqan3::seqan3)
target_link_libraries (needle robin_hood)
target_include_directories (needle PUBLIC include)

# Messages
string (ASCII 27 Esc)
set (FontBold "${Esc}[1m")
set (FontReset "${Esc}[m")

# Dependency: SeqAn3.
find_package (SeqAn3 QUIET REQUIRED HINTS lib/seqan3/build_system)

add_subdirectory (src)
message (STATUS "${FontBold}You can run `make` to build the application.${FontReset}")

## DOCUMENTATION

add_subdirectory (doc EXCLUDE_FROM_ALL)

## TEST

enable_testing ()
add_subdirectory (test EXCLUDE_FROM_ALL)
