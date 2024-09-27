<!--
SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
SPDX-License-Identifier: BSD-3-Clause
-->

# Documentation

To build the API documentation, you need to install `doxygen`. We use Doxygen version 1.8.17.
Run `make doc` and open the API documentation via `open doc/html/index.html`.

You can run `doxygen -u doxygen_cfg` to convert the configuration file to its verbose version.
`doxygen -u -s doxygen_cfg` will convert it back to its compact version.
