#pragma once
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/std/filesystem>

void convert_fastq(std::filesystem::path fastq_file, std::filesystem::path out);
