cmake_minimum_required (VERSION 3.4)
project (seqan3_tutorial CXX)
set(SeqAn3_DIR "${CMAKE_SOURCE_DIR}/../seqan3/build_system")
find_package (SeqAn3 REQUIRED)
add_executable (count
                count.cpp)
add_executable (ibf
                ibf.cpp)
add_executable (search
                search.cpp)
target_link_libraries (count seqan3::seqan3)
target_link_libraries (ibf seqan3::seqan3)
target_link_libraries (search seqan3::seqan3)
