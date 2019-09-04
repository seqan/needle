cmake_minimum_required (VERSION 3.4)
project (seqan3_tutorial CXX)
set(SeqAn3_DIR "${CMAKE_SOURCE_DIR}/../seqan3/build_system")
find_package (SeqAn3 REQUIRED)
add_executable (count_minimizer_occurences
                count_minimizer_occurences3.cpp)
add_executable (count_minimizer_occurences_transcript
                count_minimizer_occurences_transcript3.cpp)
add_executable (multiple
                multiple_count_minimizer_occurences3.cpp)
add_executable (ibf
                ibf.cpp)
add_executable (search
                query.cpp)
add_executable (testibf
                testibf.cpp)
add_executable (dist
                minimizer_occurence_dist.cpp)
target_link_libraries (count_minimizer_occurences seqan3::seqan3)
target_link_libraries (multiple seqan3::seqan3)
target_link_libraries (ibf seqan3::seqan3)
target_link_libraries (search seqan3::seqan3)
target_link_libraries (testibf seqan3::seqan3)
target_link_libraries (count_minimizer_occurences_transcript seqan3::seqan3)
target_link_libraries (dist seqan3::seqan3)
