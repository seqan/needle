# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.25...3.30)

find_package (OpenMP REQUIRED COMPONENTS CXX)
add_library ("${PROJECT_NAME}_lib" STATIC estimate.cpp ibf.cpp)
target_link_libraries ("${PROJECT_NAME}_lib" PUBLIC seqan3::seqan3)
target_link_libraries ("${PROJECT_NAME}_lib" PUBLIC robin_hood::robin_hood)
target_link_libraries ("${PROJECT_NAME}_lib" PUBLIC OpenMP::OpenMP_CXX)
target_include_directories ("${PROJECT_NAME}_lib" PUBLIC ../include)

add_executable ("${PROJECT_NAME}" main.cpp)
target_link_libraries ("${PROJECT_NAME}" PRIVATE "${PROJECT_NAME}_lib")
set_target_properties ("${PROJECT_NAME}" PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
target_compile_definitions ("${PROJECT_NAME}" PRIVATE NEEDLE_VERSION="${NEEDLE_ARGPARSE_VERSION}")
target_compile_definitions ("${PROJECT_NAME}" PRIVATE NEEDLE_DATE="${NEEDLE_ARGPARSE_DATE}")

include (GNUInstallDirs)
install (TARGETS ${PROJECT_NAME} RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")

option (NEEDLE_WITH_WERROR "Report compiler warnings as errors." ON)

if (NEEDLE_WITH_WERROR)
    target_compile_options (${PROJECT_NAME}_lib PUBLIC "-Werror")
endif ()
