// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "misc/calculate_cutoff.hpp"

#include <seqan3/io/sequence_file/format_fasta.hpp>

#include "misc/check_for_fasta_format.hpp"

// Determine cutoff for one experiment
uint8_t calculate_cutoff(std::filesystem::path sequence_file, int samples)
{
    // Cutoff according to Mantis paper -1 because we use "<" and not "<="
    uint16_t const default_cutoff{49};
    uint8_t cutoff{default_cutoff};
    std::array<uint16_t, 4> const cutoffs{0, 2, 9, 19};
    std::array<uint64_t, 4> const cutoff_bounds{314'572'800, 524'288'000, 1'073'741'824, 3'221'225'472};
    cutoff = default_cutoff;

    // Since the curoffs are based on the filesize of a gzipped fastq file, we try account for the other cases:
    // We multiply by two if we have fasta input.
    // We divide by 3 if the input is not compressed.
    bool const is_compressed = sequence_file.extension() == ".gz" || sequence_file.extension() == ".bgzf"
                            || sequence_file.extension() == ".bz2";
    bool const is_fasta = is_compressed
                            ? check_for_fasta_format(seqan3::format_fasta::file_extensions, sequence_file.stem())
                            : check_for_fasta_format(seqan3::format_fasta::file_extensions, sequence_file.extension());
    size_t const filesize =
        std::filesystem::file_size(sequence_file) * samples * (is_fasta ? 2 : 1) / (is_compressed ? 1 : 3);

    for (size_t k = 0; k < cutoff_bounds.size(); ++k)
    {
        if (filesize <= cutoff_bounds[k])
        {
            cutoff = cutoffs[k];
            break;
        }
    }
    return cutoff;
}
