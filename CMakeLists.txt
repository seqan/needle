cmake_minimum_required (VERSION 3.4)

## CUSTOMISE

# Define the application name and version.
project (my_app VERSION 1.0.0)


## BUILD

add_subdirectory(src)
add_executable ("${PROJECT_NAME}" src/main.cpp)
target_link_libraries ("${PROJECT_NAME}" "${PROJECT_NAME}_object")


## TEST

# Dependency: Googletest
add_subdirectory(lib/gtest)

enable_testing ()
add_subdirectory(test/api)
