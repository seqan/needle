// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements chopper_layout.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <chopper/chopper_layout.hpp>
#include <chopper/set_up_parser.hpp>

void chopper_layout(sharg::parser & parser)
{
    chopper::configuration config;
    set_up_parser(parser, config);
    parser.info.author = "Svenja Mehringer";
    parser.info.email = "svenja.mehringer@fu-berlin.de";

    chopper::chopper_layout(config, parser);
}
