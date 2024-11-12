# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.25...3.30)

# Define the application name and version.
project (needle
         LANGUAGES CXX
         VERSION 1.0.4
         DESCRIPTION "A fast and space-efficient pre-filter for estimating the quantification of very large collections of nucleotide sequences"
         HOMEPAGE_URL "https://github.com/seqan/needle"
)

set (NEEDLE_ARGPARSE_VERSION
     "${needle_VERSION}-rc.1"
     CACHE STRING "Needle version to display in the help message."
)
set (NEEDLE_ARGPARSE_DATE
     "2024-11-08"
     CACHE STRING "Needle's \"Last update:\" date to display in the help message."
)

# This allows including `*.cmake` files from the `cmake` directory without specifying the full path.
list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# Specify the directories where to store the built archives, libraries and executables.
include (output_directories)

# Options.
option (NEEDLE_TEST "Enable testing for needle." ON)
option (NEEDLE_DOCS "Enable documentation for needle." OFF)
option (NEEDLE_PACKAGE "Enable packaging for needle." OFF)

# Enable LTO if supported.
include (CheckIPOSupported)
check_ipo_supported (RESULT NEEDLE_HAS_LTO OUTPUT NEEDLE_HAS_LTO_OUTPUT)
if (NEEDLE_HAS_LTO)
    set (CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif ()

# Add packages.
# We use CPM for package management: https://github.com/cpm-cmake/CPM.cmake
# The added packages (e.g., hibf, sharg, seqan3) are defined in the `cmake/package-lock.cmake` file.
if (NEEDLE_PACKAGE)
    set (CPM_SOURCE_CACHE "${CMAKE_CURRENT_BINARY_DIR}/vendor")
    set (CPM_USE_LOCAL_PACKAGES OFF)
    include (package)
elseif (EXISTS "${CMAKE_CURRENT_LIST_DIR}/vendor")
    set (CPM_SOURCE_CACHE "${CMAKE_CURRENT_LIST_DIR}/vendor")
endif ()

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
    add_subdirectory (test)
endif ()

if (NEEDLE_DOCS)
    add_subdirectory (doc EXCLUDE_FROM_ALL)
endif ()
