# SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

# CPM Package Lock
# This file should be committed to version control

# cmake-format: off

# hibf
set (NEEDLE_HIBF_VERSION eca65729de77a5d411160c487adcafe8f972342b CACHE STRING "")
CPMDeclarePackage (hibf
                   NAME hibf
                   GIT_TAG ${NEEDLE_HIBF_VERSION} # main
                   GITHUB_REPOSITORY seqan/hibf
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "INSTALL_HIBF OFF"
)

# chopper
set (NEEDLE_CHOPPER_VERSION 0d03e6e80bd4f4ce82eaa8b059590de2bf7d4f7c CACHE STRING "")
CPMDeclarePackage (chopper
                   NAME chopper
                   GIT_TAG ${NEEDLE_CHOPPER_VERSION} # fast_layout_newest_only_fast_layout
                   GITHUB_REPOSITORY smehringer/chopper
                   SYSTEM TRUE
                   OPTIONS "CHOPPER_INSTALL OFF" "CHOPPER_BUILD_DOC OFF" "CHOPPER_BUILD_TEST OFF"
                           "CMAKE_MESSAGE_LOG_LEVEL WARNING"
                   EXCLUDE_FROM_ALL TRUE
)

# seqan3
set (NEEDLE_SEQAN3_VERSION 3.4.2 CACHE STRING "")
CPMDeclarePackage (seqan3
                   NAME seqan3
                   VERSION ${NEEDLE_SEQAN3_VERSION}
                   GIT_TAG ${NEEDLE_SEQAN3_VERSION}
                   GITHUB_REPOSITORY seqan/seqan3
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "INSTALL_SEQAN3 OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING"
)

# sharg
set (NEEDLE_SHARG_VERSION 1.2.2 CACHE STRING "")
CPMDeclarePackage (sharg
                   NAME sharg
                   VERSION ${NEEDLE_SHARG_VERSION}
                   GIT_TAG ${NEEDLE_SHARG_VERSION}
                   GITHUB_REPOSITORY seqan/sharg-parser
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "INSTALL_SHARG OFF" "INSTALL_TDL OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING" "BUILD_TESTING OFF"
)

# googletest
set (NEEDLE_GOOGLETEST_VERSION 1.17.0 CACHE STRING "")
CPMDeclarePackage (googletest
                   NAME GTest
                   VERSION ${NEEDLE_GOOGLETEST_VERSION}
                   GITHUB_REPOSITORY google/googletest
                   SYSTEM TRUE
                   OPTIONS "BUILD_GMOCK OFF" "INSTALL_GTEST OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING"
)

# use_ccache
set (USE_CCACHE_VERSION d2a54ef555b6fc2d496a4c9506dbeb7cf899ce37 CACHE STRING "")
CPMDeclarePackage (use_ccache
                   NAME use_ccache
                   GIT_TAG ${USE_CCACHE_VERSION} # main
                   GITHUB_REPOSITORY seqan/cmake-scripts
                   SOURCE_SUBDIR ccache
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
)

# cmake-format: on
