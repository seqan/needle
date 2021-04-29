cmake_minimum_required (VERSION 3.8)

add_library ("${PROJECT_NAME}_lib" STATIC ibf.cpp estimate.cpp)
target_link_libraries ("${PROJECT_NAME}_lib" PUBLIC seqan3::seqan3)
target_link_libraries ("${PROJECT_NAME}_lib" PUBLIC robin_hood)
find_package(OpenMP REQUIRED)
target_link_libraries("${PROJECT_NAME}_lib" PRIVATE OpenMP::OpenMP_CXX)

target_include_directories ("${PROJECT_NAME}_lib" PUBLIC ../include)

add_executable ("${PROJECT_NAME}" main.cpp)
target_link_libraries ("${PROJECT_NAME}" "${PROJECT_NAME}_lib")
