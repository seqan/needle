# Minimum cmake version
cmake_minimum_required(VERSION 3.8)

### Find doxygen and dependency to DOT tool
message (STATUS "Searching for doxygen.")
find_package (Doxygen QUIET)

if (${DOXYGEN_FOUND})

    ### Configure doc/developer targets.
    set(APP_TEMPLATE_DOXYFILE_IN ${CMAKE_SOURCE_DIR}/doc/doxygen_cfg)

    message(STATUS "Configuring doc.")

    set(APP_TEMPLATE_DOXYGEN_OUTPUT_DIR "${PROJECT_BINARY_DIR}/doc")
    set(APP_TEMPLATE_DOXYGEN_EXCLUDE_SYMBOLS "")
    set(APP_TEMPLATE_DOXYGEN_PREDEFINED_NDEBUG "")
    set(APP_TEMPLATE_DOXYGEN_ENABLED_SECTIONS "DEV")
    set(APP_TEMPLATE_DOXYGEN_EXTRACT_PRIVATE "YES")

    configure_file(${APP_TEMPLATE_DOXYFILE_IN} ${APP_TEMPLATE_DOXYGEN_OUTPUT_DIR}/Doxyfile)

    add_custom_target(doc
                        COMMAND ${DOXYGEN_EXECUTABLE}
                        WORKING_DIRECTORY ${APP_TEMPLATE_DOXYGEN_OUTPUT_DIR}
                        COMMENT "Generating (developer) API documentation with Doxygen."
                        VERBATIM)

endif ()
