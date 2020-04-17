cmake_minimum_required (VERSION 3.8)

unset(CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
unset(CMAKE_LIBRARY_OUTPUT_DIRECTORY)
unset(CMAKE_RUNTIME_OUTPUT_DIRECTORY)

# A macro that adds an api or cli test
macro (add_app_test test_filename test_alternative)

    file (RELATIVE_PATH source_file "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_LIST_DIR}/${test_filename}")

    get_filename_component (target "${source_file}" NAME_WE)

    add_executable (${target} ${test_filename})

    target_link_libraries (${target} "${PROJECT_NAME}_lib" seqan3::seqan3 gmock gmock_main gtest pthread)

    # Add main application as a dependency for command line tests
    if (test_alternative STREQUAL "CLI_TEST")
        add_dependencies(${target} "${PROJECT_NAME}")
    endif()

    get_filename_component (target_relative_path "${source_file}" DIRECTORY)

    if (target_relative_path)
        set (test_name "${target_relative_path}/${target}")
    else ()
        set (test_name "${target}")
    endif ()

    add_test (NAME "${test_name}" COMMAND ${target})

    unset (source_file)
    unset (target)
    unset (test_name)
endmacro ()

# A macro that adds an api test
macro (add_api_test test_filename)
    add_app_test(${test_filename} API_TEST)
endmacro()

# A macro that adds an cli test
macro (add_cli_test test_filename)
    add_app_test(${test_filename} CLI_TEST)
endmacro()

# Add subdirectories
add_subdirectory (api)
add_subdirectory (cli)
