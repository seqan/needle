cmake_minimum_required (VERSION 3.0)


## CUSTOMISE

# Define the application name and version.
project (my_app VERSION 1.0.0)

# The file that contains the main function.
set (file_with_main src/main.cpp)


## DEPENDENCIES

# Dependency: SeqAn3.
set (SeqAn3_DIR ${CMAKE_SOURCE_DIR}/lib/seqan3/build_system)
find_package (SeqAn3 REQUIRED)

# Dependency: Googletest
add_subdirectory(lib/gtest)


## BUILD

# An object library (without main) to be used in multiple targets.
file (GLOB productionSources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
list (REMOVE_ITEM productionSources ${CMAKE_CURRENT_SOURCE_DIR}/${file_with_main})
add_library (app_object OBJECT ${productionSources})
target_link_libraries (app_object seqan3::seqan3)
target_include_directories (app_object PRIVATE include)

# The app build target.
add_executable (${PROJECT_NAME} ${file_with_main})
target_link_libraries (${PROJECT_NAME} app_object seqan3::seqan3)
target_include_directories (${PROJECT_NAME} PRIVATE include)

# The api_test build target.
file (GLOB apiTestSources test/api/*.cpp)
add_executable (api_test ${apiTestSources})
target_link_libraries (api_test app_object seqan3::seqan3 gmock gmock_main gtest pthread)
target_include_directories (api_test PRIVATE include)


## TEST

# The test target, which executes all tests.
enable_testing ()
add_test (API_Test api_test)
