<!--
SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
SPDX-License-Identifier: BSD-3-Clause
-->

# API Test

Here are test files for API tests, i.e. the internal functions of the app are executed with different input parameters.
The test then validates whether the functions work correctly.
Attention: The default `make` target does not build tests.
Please invoke the build with `make api_test` or use `make test` to build and run all kinds of tests.
