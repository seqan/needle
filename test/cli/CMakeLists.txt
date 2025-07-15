# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.25...3.30)

add_app_test (needle_options_test.cpp)
add_app_test (ibf_options_test.cpp)
# add_app_test (insert_options_test.cpp) # Not yet supported for Needle with HIBF.
add_app_test (minimiser_options_test.cpp)
add_app_test (estimate_options_test.cpp)
add_app_test (count_options_test.cpp)
# add_app_test (delete_options_test.cpp) # Not yet supported for Needle with HIBF.
