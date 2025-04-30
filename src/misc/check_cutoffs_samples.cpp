// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "misc/check_cutoffs_samples.hpp"

#include <numeric>

// Check and set samples and cutoffs
void check_cutoffs_samples(std::vector<std::filesystem::path> const & sequence_files,
                           bool const paired,
                           std::vector<size_t> & samples,
                           std::vector<uint8_t> & cutoffs)
{
    if (paired) // If paired is true, a pair is seen as one sample
        samples.assign(sequence_files.size() / 2, 2);
    if (samples.empty()) // If no samples are given and not paired, every file is seen as one experiment
        samples.assign(sequence_files.size(), 1);
    if (cutoffs.size() == 1) // If one cutoff is given, every experiment gets this cutoff.
        cutoffs.assign(samples.size(), cutoffs[0]);

    // If sum of minimiser_args.samples is not equal to number of files, throw error
    else if (std::reduce(samples.begin(), samples.end()) != sequence_files.size())
        throw std::invalid_argument{"Incorrect command line input for multiple-samples."};
}
