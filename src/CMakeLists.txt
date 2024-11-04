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

include (GNUInstallDirs)
install (TARGETS ${PROJECT_NAME} RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")

# https://cmake.org/cmake/help/latest/module/CMakeDependentOption.html
# We usually want "-Werror" in the CI, but not when developing locally.
# Virtually all CI providers set the environment variable `CI` to some value.
# If `CI` is defined in the environment, we use `-Werror`, otherwise not.
# Caveat: If `CI` is not set, `-Dneedle_WITH_WERROR=ON` has no effect.
#         The other way around works. If `CI` is set, `-Dneedle_WITH_WERROR=OFF` is honored.
# Solutions: * `export CI=1` (any value after `=` works, the variable only needs to be not empty)
#            * `-DCMAKE_CXX_FLAGS="-Werror"`
include (CMakeDependentOption)
cmake_dependent_option (${PROJECT_NAME}_WITH_WERROR "Report compiler warnings as errors." ON "DEFINED ENV{CI}" OFF)

if (${PROJECT_NAME}_WITH_WERROR)
    target_compile_options (${PROJECT_NAME}_lib PUBLIC "-Werror")
endif ()
