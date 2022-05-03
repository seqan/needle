// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/needle/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/core/debug_stream.hpp>

#include "shared.h"
#include "ibf.h"
#include "estimate.h"

uint32_t w_size;
uint64_t shape{};
uint64_t se;

void initialise_min_arguments(seqan3::argument_parser & parser, min_arguments & args)
{
    parser.add_option(args.k, 'k', "kmer", "Define k-mer size for the minimisers. Default: 20.");
    parser.add_option(w_size, 'w', "window", "Define window size for the minimisers. Default: 60.");
    parser.add_option(shape, '\0', "shape", "Define a shape for the minimisers by the decimal of a bitvector, where 0 symbolizes a "
                                           "position to be ignored, 1 a position considered. Default: ungapped.");
    parser.add_option(se, '\0', "seed", "Define seed for the minimisers.");
    parser.add_option(args.path_out, 'o', "out", "Directory, where output files should be saved.");
    parser.add_option(args.threads, 't', "threads", "Number of threads to use. Default: 1.");
}

void initialise_arguments_ibf(seqan3::argument_parser & parser, estimate_ibf_arguments & ibf_args, size_t & num_hash,
                              std::vector<double> & fpr)
{
    parser.add_flag(ibf_args.compressed, 'c', "compressed", "If c is set, the IBFS are compressed. Default: Not compressed.");
    parser.add_option(fpr, 'f', "fpr", "List of bin false positive rate per expression level. If only one is given"
                                                          ", then that fpr is used for all expression levels.");
    parser.add_option(ibf_args.expression_thresholds, 'e', "expression_thresholds", "Which expression thresholds should be used for"
                                                                                " constructing the IBFs.");
    parser.add_option(ibf_args.number_expression_thresholds, 'l', "number_expression_thresholds", "Number of expression thresholds. "
                                                              "Can be set alternatively to expression_thresholds, then "
                                                              "the expression thresholds are determined automatically.");
    parser.add_option(num_hash, 'n', "hash", "Number of hash functions that should be used when constructing "
                                             "one IBF.");
}

void parsing(seqan3::argument_parser & parser, min_arguments & args)
{
    w_size = args.w_size.get();
    se = args.s.get();
    parser.parse();
    args.w_size = seqan3::window_size{w_size};
    if (shape == 0)
            args.shape = seqan3::ungapped{args.k};
    else
            args.shape = seqan3::bin_literal{shape};
    args.s = seqan3::seed{adjust_seed(args.k, se)};
}

// Initialize arguments for ibf and minimiser
void initialise_arguments_minimiser(seqan3::argument_parser & parser, minimiser_arguments & minimiser_args, std::vector<uint8_t> & cutoffs)
{
    parser.add_option(minimiser_args.include_file, '\0', "include", "Sequence file containing minimizers, only those "
                                                                    "minimizers will be considered.");
    parser.add_option(minimiser_args.exclude_file, '\0', "exclude", "Sequence file containing minimizers that should "
                                                                    "not be stored.");
    parser.add_option(minimiser_args.samples, '\0', "samples", "Define which samples belong together, sum has to be "
                                                               "equal to number of sequence files. Default: Every"
                                                               " sequence file is one sample from one experiment.");
    parser.add_flag(minimiser_args.paired, 'p', "paired", "If set, experiments are paired. Default: Not paired.");
    parser.add_option(cutoffs, '\0', "cutoff", "Define for each sample, what number of found minimisers "
                                                              "should be considered the result of a sequencing error "
                                                              "and therefore be ignored. Default: Every sample has an"
                                                              "automatically genereated cutoff, which is based on the "
                                                              "file size.");

}

void read_input_file_list(std::vector<std::filesystem::path> & sequence_files, std::filesystem::path input_file)
{
    std::ifstream fin{input_file};

    if (!fin.good() || !fin.is_open())
        throw std::runtime_error{"Could not open file " + input_file.string() + " for reading."};

    std::string line;
    while (std::getline(fin, line))
    {
        sequence_files.push_back(line);
    }
}

int run_needle_count(seqan3::argument_parser & parser)
{
    min_arguments args;
    initialise_min_arguments(parser, args);
    std::vector<std::filesystem::path> sequence_files{};
    std::filesystem::path genome_file;
    std::filesystem::path exclude_file{""};
    std::filesystem::path out_path = "./";
    bool paired = false;

    parser.info.short_description = "Get expression value depending on minimizers. This function is only used to test "
                                    "the validity of Needle's estimation approach. It estimates the expression value "
                                    "for all sequences in the genome file based on the exact minimiser occurrences of "
                                    "the given sequence files.";
    parser.add_positional_option(sequence_files, "Please provide at least one sequence file.");
    parser.add_option(genome_file, 'g', "genome", "Please provide one sequence file with transcripts.");
    parser.add_option(exclude_file, '\0', "exclude", "Please provide one sequence file with minimizers to ignore.");
    parser.add_flag(paired, 'p', "paired", "If set, experiments are paired. Default: Not paired.");

    try
    {
        parsing(parser, args);
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        seqan3::debug_stream << "Error. Incorrect command line input for count. " << ext.what() << "\n";
        return -1;
    }

    count(args, sequence_files, genome_file, exclude_file, paired);

    return 0;
}

int run_needle_estimate(seqan3::argument_parser & parser)
{
    estimate_ibf_arguments args{};
    estimate_arguments estimate_args{};
    parser.info.short_description = "Estimate expression value of transcript based on the Needle index.";
    parser.info.version = "1.0.0";
    parser.info.author = "Mitra Darvish";

    args.path_out = "expressions.out";

    parser.add_positional_option(estimate_args.search_file, "Please provide a sequence file.");
    parser.add_option(estimate_args.path_in, 'i', "in", "Directory where input files can be found.");
    parser.add_option(args.path_out, 'o', "out", "Directory, where output files should be saved.");
    parser.add_flag(estimate_args.normalization_method, 'm', "normalization-mode",
                                                            "Set, if normalization is wanted. Normalization is achieved by"
                                                            "dividing the expression value with the expression threshold of the first"
                                                            " ibf. Only make sense if every bin has its own expression "
                                                            "thresholds (which is the case if expression thresholds "
                                                            "were generated automatically)."
                                                            "Default: False.");

    try
    {
        parsing(parser, args);
    }
    catch (seqan3::argument_parser_error const & ext)                     // catch user errors
    {
        seqan3::debug_stream << "Error. Incorrect command line input for estimate. " << ext.what() << "\n";
        return -1;
    }

    call_estimate(args, estimate_args);

    return 0;
}


int run_needle_ibf(seqan3::argument_parser & parser)
{
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    std::filesystem::path input_file{};
    std::vector<std::filesystem::path> sequence_files{};
    size_t num_hash{1}; // Number of hash functions to use, default 1
    std::filesystem::path expression_by_genome_file = "";
    std::vector<double> fpr{}; // The fpr of one IBF, can be different for different expression levels
    std::vector<uint8_t> cutoffs{};

    initialise_min_arguments(parser, ibf_args);
    initialise_arguments_ibf(parser, ibf_args, num_hash, fpr);
    initialise_arguments_minimiser(parser, minimiser_args, cutoffs);

    parser.info.short_description = "Constructs the Needle index.";

    parser.add_option(input_file, 'i', "in", "Please provide one file containing a list of sequence files.");
    parser.add_positional_option(sequence_files, "Please provide at least one sequence file OR provide one file "
                                                 "containing all sequence files with -i.");
    parser.add_option(minimiser_args.experiment_names, '\0', "experiment-names", "If set, names of the experiments are stored"
                                                                                 " in a txt file.");
    parser.add_option(expression_by_genome_file, '\0', "levels-by-genome", "Sequence file containing minimizers, only "
                                                                            "those minimizers will be considered for "
                                                                            "determining the expression thresholds.");

    try
    {
        parsing(parser, ibf_args);
        if (input_file != "")
            read_input_file_list(sequence_files, input_file);
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        seqan3::debug_stream << "Error. Incorrect command line input for ibf. " << ext.what() << "\n";
        return -1;
    }

    try
    {
        ibf(sequence_files, ibf_args, minimiser_args, fpr, cutoffs, expression_by_genome_file, num_hash);
    }
    catch (const std::invalid_argument & e)
    {
        std::cerr << e.what() << std::endl;
        return -1;
    }

    return 0;
}

int run_needle_ibf_min(seqan3::argument_parser & parser)
{
    estimate_ibf_arguments ibf_args{};
    std::vector<std::filesystem::path> minimiser_files{};
    size_t num_hash{1}; // Number of hash functions to use, default 1
    std::filesystem::path expression_by_genome_file = "";
    std::vector<double> fpr{}; // The fpr of one IBF, can be different for different expression levels
    std::filesystem::path input_file{};

    parser.info.short_description = "Constructs the Needle index from the minimiser files created by needle minimiser.";

    parser.add_option(input_file, 'i', "in", "Please provide one file containing a list of sequence files.");
    parser.add_positional_option(minimiser_files, "Please provide at least one minimiser file OR provide one file "
                                                 "containing all minimiser files with -i.");

    parser.add_option(ibf_args.path_out, 'o', "out", "Directory, where output files should be saved.");
    parser.add_option(ibf_args.threads, 't', "threads", "Number of threads to use. Default: 1.");
    parser.add_option(expression_by_genome_file, '\0', "levels-by-genome", "Sequence file containing minimizers, only "
                                                                            "those minimizers will be considered for "
                                                                            "determining the expression thresholds.");

    initialise_arguments_ibf(parser, ibf_args, num_hash, fpr);

    try
    {
        parsing(parser, ibf_args);
        if (input_file != "")
            read_input_file_list(minimiser_files, input_file);
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        seqan3::debug_stream << "Error. Incorrect command line input for ibfmin. " << ext.what() << "\n";
        return -1;
    }

    try
    {
        ibf(minimiser_files, ibf_args, fpr, expression_by_genome_file, num_hash);
    }
    catch (const std::invalid_argument & e)
    {
        std::cerr << e.what() << std::endl;
        return -1;
    }

    return 0;
}

int run_needle_minimiser(seqan3::argument_parser & parser)
{
    min_arguments args{};
    minimiser_arguments minimiser_args{};
    std::vector<std::filesystem::path> sequence_files{};
    std::vector<uint8_t> cutoffs{};
    initialise_min_arguments(parser, args);
    initialise_arguments_minimiser(parser, minimiser_args, cutoffs);
    std::filesystem::path input_file{};

    parser.info.short_description = "Calculates minimiser for given experiments.";
    parser.add_option(input_file, 'i', "in", "Please provide one file containing a list of sequence files.");
    parser.add_positional_option(sequence_files, "Please provide at least one sequence file OR provide one file "
                                                 "containing all sequence files with -i.");

    try
    {
        parsing(parser, args);
        if (input_file != "")
            read_input_file_list(sequence_files, input_file);
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        seqan3::debug_stream << "Error. Incorrect command line input for minimiser. " << ext.what() << "\n";
        return -1;
    }
    try
    {
        minimiser(sequence_files, args, minimiser_args, cutoffs);
    }
    catch (const std::invalid_argument & e)
    {
        std::cerr << e.what() << std::endl;
        return -1;
    }

    return 0;
}

int main(int argc, char const ** argv)
{
    seqan3::argument_parser needle_parser{"needle", argc, argv, seqan3::update_notifications::on,
    {"count", "estimate", "ibf", "ibfmin", "minimiser"}};
    needle_parser.info.description.push_back("Needle allows you to build an Interleaved Bloom Filter (IBF) with the "
                                             "command ibf or estimate the expression of transcripts with the command "
                                             "estimate.");
    needle_parser.info.version = "1.0.0";
    needle_parser.info.author = "Mitra Darvish";

    try
    {
        needle_parser.parse(); // trigger command line parsing
    }
    catch (seqan3::argument_parser_error const & ext) // catch user errors
    {
        seqan3::debug_stream << "Error. Incorrect command. See needle help for more information." << ext.what() << "\n";
        return -1;
    }
    seqan3::argument_parser & sub_parser = needle_parser.get_sub_parser(); // hold a reference to the sub_parser
    if (sub_parser.info.app_name == std::string_view{"needle-count"})
        run_needle_count(sub_parser);
    else if (sub_parser.info.app_name == std::string_view{"needle-estimate"})
        run_needle_estimate(sub_parser);
    else if (sub_parser.info.app_name == std::string_view{"needle-ibf"})
        run_needle_ibf(sub_parser);
    else if (sub_parser.info.app_name == std::string_view{"needle-ibfmin"})
        run_needle_ibf_min(sub_parser);
    else if (sub_parser.info.app_name == std::string_view{"needle-minimiser"})
        run_needle_minimiser(sub_parser);
}
