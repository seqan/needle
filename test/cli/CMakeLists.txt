cmake_minimum_required (VERSION 3.8)

add_cli_test (fastq_to_fasta_options_test.cpp)
target_use_datasources (fastq_to_fasta_options_test FILES in.fastq)
