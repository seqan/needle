cmake_minimum_required (VERSION 3.4)

## DEPENDENCIES

# Dependency: SeqAn3.
find_package (SeqAn3 REQUIRED HINTS ../lib/seqan3/build_system)


## BUILD

# An object library (without main) to be used in multiple targets.
add_library ("${PROJECT_NAME}_object" OBJECT give_me_five.cpp)
target_link_libraries ("${PROJECT_NAME}_object" PUBLIC seqan3::seqan3)
target_include_directories ("${PROJECT_NAME}_object" PUBLIC ../include)
