cmake_minimum_required (VERSION 3.8)

unset(CMAKE_RUNTIME_OUTPUT_DIRECTORY)

# Get a specific component of a test file which follows the seqan3 naming scheme.
# e.g. target_source_file = "range/views/take.cpp"
# component:
#  * TARGET_NAME - the target name (e.g. "take")
#  * TEST_NAME - the test name which includes the target_path (e.g. "range/views/take")
macro (seqan3_test_component VAR target_source_file component_name_)
    string (TOUPPER "${component_name_}" component_name)

    get_filename_component (target_relative_path "${target_source_file}" DIRECTORY)
    get_filename_component (target_name "${target_source_file}" NAME_WE)

    if (component_name STREQUAL "TARGET_NAME")
        set (${VAR} "${target_name}")
    elseif (component_name MATCHES "TEST_NAME")
        if (target_relative_path)
            set (${VAR} "${target_relative_path}/${target_name}")
        else ()
            set (${VAR} "${target_name}")
        endif ()
    endif ()

    unset (target_name)
    unset (target_relative_path)
endmacro ()

# A macro that adds an api or cli test
macro (add_app_test unit_test_cpp)
    file (RELATIVE_PATH unit_test "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_LIST_DIR}/${unit_test_cpp}")
    seqan3_test_component (target "${unit_test}" TARGET_NAME)
    seqan3_test_component (test_name "${unit_test}" TEST_NAME)

    add_executable (${target} ${unit_test_cpp})

    target_link_libraries (${target} "${PROJECT_NAME}_lib" seqan3::seqan3 gmock gmock_main gtest pthread)
    add_dependencies(${target} "${PROJECT_NAME}")
    add_test (NAME "${test_name}" COMMAND ${target})

    unset (unit_test)
    unset (target)
    unset (test_name)
endmacro ()

# Add subdirectories
add_subdirectory (api)
add_subdirectory (cli)

