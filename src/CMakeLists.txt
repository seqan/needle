cmake_minimum_required (VERSION 3.25...3.30)

find_package (OpenMP REQUIRED COMPONENTS CXX)
add_library ("${PROJECT_NAME}_lib" STATIC estimate.cpp ibf.cpp)
target_link_libraries ("${PROJECT_NAME}_lib" PUBLIC seqan3::seqan3)
target_link_libraries ("${PROJECT_NAME}_lib" PUBLIC robin_hood::robin_hood)
target_link_libraries ("${PROJECT_NAME}_lib" PUBLIC OpenMP::OpenMP_CXX)
target_include_directories ("${PROJECT_NAME}_lib" PUBLIC ../include)

add_executable ("${PROJECT_NAME}" main.cpp)
target_link_libraries ("${PROJECT_NAME}" PRIVATE "${PROJECT_NAME}_lib")
