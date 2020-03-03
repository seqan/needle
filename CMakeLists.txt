cmake_minimum_required (VERSION 3.4)

## CUSTOMISE

# Define the application name and version.
project (needle VERSION 1.0.0)

## BUILD

add_subdirectory(src)
set(DATA_DIR DATADIR=${CMAKE_CURRENT_SOURCE_DIR})
add_compile_definitions(DATA_DIR="${CMAKE_CURRENT_LIST_DIR}/test/data")
set(CMAKE_BUILD_TYPE Release)
add_executable (needle src/needle.cpp)
target_link_libraries (needle seqan3::seqan3)
target_include_directories (needle PUBLIC include)

## TEST

# Dependency: Googletest

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


# Now simply link against gtest or gtest_main as needed. Eg
#target_link_libraries(test-needle gtest_main)
# add_test(NAME example_test COMMAND example)
enable_testing()
add_subdirectory(test/api)
