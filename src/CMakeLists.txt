cmake_minimum_required (VERSION 3.9)

find_package(OpenMP REQUIRED)
add_library ("${PROJECT_NAME}_lib" STATIC ibf.cpp estimate.cpp)
target_link_libraries ("${PROJECT_NAME}_lib" PUBLIC
    seqan3::seqan3
    robin_hood
    OpenMP::OpenMP_CXX
    sharg::sharg
)

target_include_directories ("${PROJECT_NAME}_lib" PUBLIC ../include)

add_executable ("${PROJECT_NAME}" main.cpp)
target_link_libraries ("${PROJECT_NAME}" PRIVATE "${PROJECT_NAME}_lib")
