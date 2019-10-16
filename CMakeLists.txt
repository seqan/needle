cmake_minimum_required (VERSION 3.4)
project (seqan3_tutorial CXX)
set(SeqAn3_DIR "${SEQAN_DIR}/build_system")
find_package (SeqAn3 REQUIRED)
add_executable (count
                count.cpp)
add_executable (needle-ibf
                ibf.cpp)
add_executable (needle-search
                search.cpp)
target_link_libraries (count seqan3::seqan3)
target_link_libraries (needle-ibf seqan3::seqan3)
target_link_libraries (needle-search seqan3::seqan3)
