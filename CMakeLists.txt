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

# An option to disable configuring and building the tests. Tests are enabled by default.
# If your project-name (line 8 of this file) is `app-template`, the option will be `app-template_TEST`.
# It can be used when calling CMake: `cmake .. -Dapp-template_TEST=OFF`.
# It is good practice to allow disabling tests. If another project includes your application,
# it might not want to build your tests.
option (${PROJECT_NAME}_TEST "Enable testing for ${PROJECT_NAME}." ON)

if (${PROJECT_NAME}_TEST)
    # Add the tests. This will include `test/CMakeLists.txt`.
    add_subdirectory (test EXCLUDE_FROM_ALL)
endif ()

# add_subdirectory (doc EXCLUDE_FROM_ALL)

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

# set (CPACK_GENERATOR "TXZ")

# set (CPACK_PACKAGE_VERSION "${needle_VERSION}")
# set (CPACK_PACKAGE_VENDOR "seqan")
# set (CPACK_PACKAGE_CHECKSUM "SHA256")
# set (CPACK_RESOURCE_FILE_LICENSE "${needle_SOURCE_DIR}/LICENSE.md")
# set (CPACK_RESOURCE_FILE_README "${needle_SOURCE_DIR}/README.md")

# set (CPACK_SOURCE_GENERATOR "TXZ")
# set (CPACK_SOURCE_IGNORE_FILES "\\\\.git($|/)")

# include (CPack)
