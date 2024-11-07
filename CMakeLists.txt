# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.25...3.30)

# Define the application name and version.
project (needle
         LANGUAGES CXX
         VERSION 1.0.3
         DESCRIPTION "A fast and space-efficient pre-filter for estimating the quantification of very large collections of nucleotide sequences"
)

# This allows including `*.cmake` files from the `cmake` directory without specifying the full path.
list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# Specify the directories where to store the built archives, libraries and executables.
include (output_directories)

# An option to disable configuring and building the tests. Tests are enabled by default.
# If your project-name (line 8 of this file) is `app-template`, the option will be `app-template_TEST`.
# It can be used when calling CMake: `cmake .. -Dapp-template_TEST=OFF`.
# It is good practice to allow disabling tests. If another project includes your application,
# it might not want to build your tests.
option (NEEDLE_TEST "Enable testing for NEEDLE." ON)
option (NEEDLE_DOCS "Enable documentation for NEEDLE." OFF)
option (NEEDLE_PACKAGE "Enable packaging for NEEDLE." OFF)

if (NEEDLE_PACKAGE)
    set (CPM_SOURCE_CACHE "${CMAKE_CURRENT_BINARY_DIR}/vendor")
    set (CPM_USE_LOCAL_PACKAGES OFF)
    include (package)
elseif (EXISTS "${CMAKE_CURRENT_LIST_DIR}/vendor")
    set (CPM_SOURCE_CACHE "${CMAKE_CURRENT_LIST_DIR}/vendor")
endif ()

# Add packages.
# We use CPM for package management: https://github.com/cpm-cmake/CPM.cmake
# The added packages (e.g., hibf, sharg, seqan3) are defined in the `cmake/package-lock.cmake` file.
include (CPM)
CPMUsePackageLock (${CMAKE_CURRENT_LIST_DIR}/cmake/package-lock.cmake)

# Use ccache if available. This speeds up the build process by caching files that have been compiled before.
CPMGetPackage (use_ccache)
CPMGetPackage (seqan3)
CPMGetPackage (robin-hood)

# Add the application. This will include `src/CMakeLists.txt`.
add_subdirectory (src)

if (NEEDLE_TEST)
    list (APPEND CMAKE_CTEST_ARGUMENTS "--output-on-failure;--no-tests=error") # Must be before `enable_testing ()`.
    enable_testing ()
    add_subdirectory (test EXCLUDE_FROM_ALL)
endif ()

if (NEEDLE_DOCS)
    add_subdirectory (doc EXCLUDE_FROM_ALL)
endif ()
