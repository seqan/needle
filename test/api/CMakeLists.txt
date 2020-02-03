cmake_minimum_required (VERSION 3.4)

# The api_test build target.
add_executable (give_me_five_test give_me_five_test.cpp)
target_link_libraries (give_me_five_test "${PROJECT_NAME}_object" gmock gmock_main gtest pthread)
add_test (give_me_five_test give_me_five_test)
