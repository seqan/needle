# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

cmake_minimum_required (VERSION 3.25)

# This includes `cmake/test/config.cmake` which takes care of setting up the test infrastructure. It also provides
# the `add_app_test` macro, which is used to add tests to the CMake test suite.
include (test/config)

# This includes `test/data/datasources.cmake`, which makes test data available to the tests.
include (data/datasources.cmake)

add_subdirectory (api)
add_subdirectory (cli)

message (STATUS "You can run `make tests` to build tests.")
message (STATUS "You can run `make check` to build and run tests.")
