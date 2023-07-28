cmake_minimum_required (VERSION 3.9)

# The api_test build target.
add_api_test (count_test.cpp)
add_api_test (insert_delete_test.cpp)
add_api_test (estimate_test.cpp)
add_api_test (ibf_test.cpp)
add_api_test (ibfmin_test.cpp)
add_api_test (minimiser_test.cpp)
