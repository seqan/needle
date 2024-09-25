# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

set (DATASOURCES_DATA_DIR "${${PROJECT_NAME}_SOURCE_DIR}/test/data")

file (GLOB_RECURSE datasources
      LIST_DIRECTORIES false
      RELATIVE ${DATASOURCES_DATA_DIR}
      CONFIGURE_DEPENDS ${DATASOURCES_DATA_DIR}/*
)
list (REMOVE_ITEM datasources datasources.cmake README.md)
list (FILTER datasources EXCLUDE REGEX "\.license")

foreach (datasource IN LISTS datasources)
    get_filename_component (datasource_name "${datasource}" NAME)
    file (SHA256 ${DATASOURCES_DATA_DIR}/${datasource} datasource_hash)

    declare_datasource (FILE ${datasource_name}
                        URL ${DATASOURCES_DATA_DIR}/${datasource}
                        URL_HASH SHA256=${datasource_hash}
    )
endforeach ()
