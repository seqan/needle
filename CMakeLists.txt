cmake_minimum_required (VERSION 3.4)

# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                 EXCLUDE_FROM_ALL)

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
  include_directories("${gtest_SOURCE_DIR}/include")
endif()

set (EXECUTABLE_OUTPUT_PATH ${OUTPUT_PATH})
set (LIBRARY_OUTPUT_PATH ${OUTPUT_PATH})
project (seqan3_tutorial CXX)
set(SeqAn3_DIR "${SEQAN_DIR}/build_system")
set(CMAKE_BUILD_TYPE Release)
find_package (SeqAn3 REQUIRED)
add_executable (count
                count.cpp)
add_executable (ibf
                needle.cpp)
add_executable (needle-search
                search.cpp)
add_executable (test
                test.cpp)
target_link_libraries (count seqan3::seqan3)
target_link_libraries (ibf seqan3::seqan3)
target_link_libraries (needle-search seqan3::seqan3)
target_link_libraries (test seqan3::seqan3)
# Now simply link against gtest or gtest_main as needed. Eg
target_link_libraries(test gtest_main)
# add_test(NAME example_test COMMAND example)
