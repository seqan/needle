cmake_minimum_required (VERSION 3.4)

## DEPENDENCIES

# Dependency: SeqAn3.
find_package (SeqAn3 REQUIRED HINTS ../lib/seqan3/build_system)
add_executable (count count.cpp)
target_link_libraries (count seqan3::seqan3)
target_include_directories (count PUBLIC ../include)
