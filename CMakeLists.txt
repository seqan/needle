cmake_minimum_required (VERSION 3.4)
project (seqan3_tutorial CXX)
set(SeqAn3_DIR "${SEQAN_DIR}/build_system")
set(CMAKE_BUILD_TYPE Release)
find_package (SeqAn3 REQUIRED)
add_executable (count-mini
                count_minimizer.cpp)
target_link_libraries (count-mini seqan3::seqan3)
