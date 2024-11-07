# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

## PACKAGE
# Change in top-level CMakeLists.txt:
#   Version -> project (needle VERSION x.x.x)
#   Date -> set (NEEDLE_ARGPARSE_DATE "YYYY-MM-DD")
# No dependencies should be locally installed.
# To package, create a clean directory:
#
# mkdir needle-package
# cd needle-package
# git clone https://github.com/seqan/needle
# cd needle
# cd ../../../
# mkdir package
# cd package
# cmake ../needle
# cmake --build . --target package_source
#
# Will create needle-[VERSION]-Source.tar.xz{,.sha256}.
#
# Alternatively, the ci_install workflow will export the source package as artifact.

cmake_minimum_required (VERSION 3.7...3.30)

set (CPACK_GENERATOR "TXZ")

set (CPACK_PACKAGE_VERSION "${needle_VERSION}")
set (CPACK_PACKAGE_VENDOR "seqan")
set (CPACK_PACKAGE_CHECKSUM "SHA256")
set (CPACK_RESOURCE_FILE_LICENSE "${needle_SOURCE_DIR}/LICENSE.md")
set (CPACK_RESOURCE_FILE_README "${needle_SOURCE_DIR}/README.md")

# Already being called on source package, i.e. CPM is already downloaded.
if (NOT CPM_DOWNLOAD_LOCATION)
    set (CPM_DOWNLOAD_LOCATION "${needle_SOURCE_DIR}/cmake/CPM.cmake")
endif ()

configure_file ("${needle_SOURCE_DIR}/cmake/cpack_install.cmake.in" "${CMAKE_CURRENT_BINARY_DIR}/cpack_install.cmake"
                @ONLY
)
set (CPACK_INSTALL_SCRIPT "${CMAKE_CURRENT_BINARY_DIR}/cpack_install.cmake")

# Source Package
set (CPACK_SOURCE_GENERATOR "TXZ")
list (APPEND CPACK_SOURCE_IGNORE_FILES "/\.git($|/)")
list (APPEND CPACK_SOURCE_IGNORE_FILES "/\.github/")
list (APPEND CPACK_SOURCE_IGNORE_FILES "/\.vscode/")
list (APPEND CPACK_SOURCE_IGNORE_FILES "/build/")

include (CPack)
