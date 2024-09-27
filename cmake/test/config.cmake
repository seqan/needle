# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

include (test/coverage)

CPMGetPackage (googletest)

list (APPEND CMAKE_CTEST_ARGUMENTS "--output-on-failure") # Must be before `enable_testing ()`.
enable_testing ()

# Set directories for test output files, input data and binaries.
file (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output)
add_definitions (-DOUTPUTDIR=\"${CMAKE_CURRENT_BINARY_DIR}/output/\")
add_definitions (-DDATADIR=\"${CMAKE_CURRENT_BINARY_DIR}/data/\")
add_definitions (-DBINDIR=\"${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/\")
add_definitions (-DAPPNAME=\"${PROJECT_NAME}\")

# Add the test interface library.
if (NOT TARGET ${PROJECT_NAME}_test)
    add_library (${PROJECT_NAME}_test INTERFACE)
    target_compile_options (${PROJECT_NAME}_lib PUBLIC "-pedantic" "-Wall" "-Wextra" "-Werror")

    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        # Disable bogus warnings in GCC12.
        if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 12 AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 13)
            target_compile_options (${PROJECT_NAME}_test INTERFACE "-Wno-array-bounds" "-Wno-stringop-overread")
        endif ()

        # Warn about failed return value optimization.
        if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 14)
            target_compile_options (${PROJECT_NAME}_lib PUBLIC "-Wnrvo")
        endif ()
    endif ()

    target_link_libraries (${PROJECT_NAME}_test INTERFACE "${PROJECT_NAME}_lib" "GTest::gtest_main")

    # !Workaround: Get seqan3 test include dir from seqan3 target
    find_path (SEQAN3_TEST_INCLUDE_DIR
               NAMES seqan3/test/tmp_directory.hpp
               HINTS "${seqan3_SOURCE_DIR}/test/include"
    )
    target_include_directories (${PROJECT_NAME}_test SYSTEM INTERFACE "${SEQAN3_TEST_INCLUDE_DIR}")

    add_library (${PROJECT_NAME}::test ALIAS ${PROJECT_NAME}_test)
endif ()

# Add the check target that builds and runs tests.
add_custom_target (check COMMAND ${CMAKE_CTEST_COMMAND} ${CMAKE_CTEST_ARGUMENTS})
add_custom_target (tests)

macro (add_app_test test_filename)
    file (RELATIVE_PATH source_file "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_LIST_DIR}/${test_filename}")
    get_filename_component (target "${source_file}" NAME_WE)

    add_executable (${target} ${test_filename})
    target_link_libraries (${target} ${PROJECT_NAME}::test)

    add_dependencies (${target} ${PROJECT_NAME})
    add_dependencies (check ${target})
    add_dependencies (tests ${target})

    add_test (NAME ${target} COMMAND ${target})

    unset (source_file)
    unset (target)
endmacro ()
