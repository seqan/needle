// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "misc/debug.hpp"

#include <iostream>

void log(std::string_view const message, std::source_location const location)
{
    std::clog << location.file_name() << ':' << location.line() << ':' << location.column() << " `"
              << location.function_name() << "`: " << message << '\n';
}
