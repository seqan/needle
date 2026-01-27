# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.25...3.30)

if (TARGET "${PROJECT_NAME}_layout")
    return ()
endif ()

add_library ("${PROJECT_NAME}_layout" STATIC layout.cpp)
target_link_libraries ("${PROJECT_NAME}_layout" PRIVATE "chopper::chopper")
