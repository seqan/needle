#pragma once
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/std/filesystem>

/*! \brief Function, converting fastq files to fasta files.
 *  \param fastq_file input file path to the fastq file
 *  \param out output file path for the fasta file
 *
 *  Simple function, converting fastq files to fasta files using the seqan3 library.
 *  For more information about the SeqAn Library functions see https://docs.seqan.de/seqan/3-master-user/.
 */
void convert_fastq(std::filesystem::path fastq_file, std::filesystem::path out);
