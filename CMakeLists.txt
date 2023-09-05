cmake_minimum_required (VERSION 3.10)

## CUSTOMISE

# Define the application name and version.
project (needle VERSION 1.0.3)

## BUILD

# Make Release default build type
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release CACHE STRING
         "Choose the type of build, options are: Debug Release Coverage RelWithDebInfo MinSizeRel."
         FORCE)
endif ()

# Specify the directories where to store the built archives, libraries and executables
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Messages
string (ASCII 27 Esc)
set (FontBold "${Esc}[1m")
set (FontReset "${Esc}[m")

# Dependency: SeqAn3.
find_package (SeqAn3 QUIET REQUIRED HINTS lib/seqan3/build_system)

# Use ccache.
include ("${SEQAN3_CLONE_DIR}/test/cmake/seqan3_require_ccache.cmake")
seqan3_require_ccache ()

# Dependency: robin_hood.
add_subdirectory(lib/robin-hood-hashing)

# Add the application.
add_subdirectory (src)
message (STATUS "${FontBold}You can run `make` to build the application.${FontReset}")


## DOCUMENTATION

add_subdirectory (doc EXCLUDE_FROM_ALL)

## TEST

enable_testing ()
add_subdirectory (test EXCLUDE_FROM_ALL)

## PACKAGE
# Change version in this file: project (needle VERSION x.x.x)
# Change version in src/main.cpp: parser.info.version = "x.x.x";
# To package, create a clean directory:
#
# mkdir needle-package
# cd needle-package
# git clone https://github.com/seqan/needle
# cd needle
# # git checkout VERSION # optional: Checkout a tag if it exists.
# git submodule update --init
# cd lib/seqan3
# git submodule update --init
# cd ../../../
# mkdir package
# cd package
# cmake ../needle
# cmake --build . --target package_source
#
# Will create needle-[VERSION]-Source.tar.xz{,.sha256}.
# We do git submodule update --init instead of a recursive clone, because we do not want to package the recursive
# submodules, but only what we need. The script handles all submodules at the time of writing, but newer submodules
# might be missing.
# This should be done in a clean directory, because everything in the needle checkout will be packaged, i.e., also
# local build directories if they exist.

set (CPACK_GENERATOR "TXZ")

set (CPACK_PACKAGE_VERSION "${needle_VERSION}")
set (CPACK_PACKAGE_VENDOR "seqan")
set (CPACK_PACKAGE_CHECKSUM "SHA256")
set (CPACK_RESOURCE_FILE_LICENSE "${needle_SOURCE_DIR}/LICENSE.md")
set (CPACK_RESOURCE_FILE_README "${needle_SOURCE_DIR}/README.md")

set (CPACK_SOURCE_GENERATOR "TXZ")
set (CPACK_SOURCE_IGNORE_FILES "\\\\.git($|/)")

include (CPack)
