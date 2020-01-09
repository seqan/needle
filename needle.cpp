#include <seqan3/argument_parser/all.hpp>
#include "ibf.h"
#include "search.h"

// =====================================================================================================================
// ibf
// =====================================================================================================================
/*
int main(int argc, char const ** argv) //run_needle_ibf(argument_parser & parser)
{
    seqan3::argument_parser parser("needle-ibf", argc, argv);
    cmd_arguments args{};
    parser.info.author = "Mitra Darvish";
    parser.info.short_description = "Constructs IBF on Minimizer.";
    parser.add_positional_option(args.sequence_files, "Please provide at least one sequence file.");
    parser.add_option(args.genome_file, 'g', "genom-mask", "Genom file used as a mask.");
    parser.add_option(args.bits, 'z', "size", "List of sizes in bits for IBF per expression rate.");
    parser.add_option(args.num_hash, 'n', "hash", "Number of hash functions that should be used when constructing "
                      "one IBF.");
    parser.add_option(args.path_out, 'o', "out", "Directory, where output files should be saved.");
    parser.add_option(args.expression_levels, 'e', "expression_levels", "Which expression levels should be used for "
                      "constructing the IBFs.");
    parser.add_option(args.samples, 'm', "multiple-samples", "Define which samples belong together, sum has to be equal"
                      " to number of sequence files. Default: Every sequence file is one sample from one experiment.");
    parser.add_option(args.aggregate_by, 'a', "aggregate-by", "Choose your method of aggregation: mean, median or "
                      "random. Default: median.");
    parser.add_option(args.random, 'r', "random-samples", "Choose the number of random sequences to pick from when "
                      "using aggregation method random. Default: 1000.");

    parser.add_flag(args.compressed, 'c', "compressed", "If set ibf is compressed. Default: Not compressed.");
    parser.add_option(args.k, 'k', "kmer", "Define kmer size.");
    parser.add_option(args.window_size, 'w', "window", "Define window size.");
    parser.add_option(args.shape, 'p', "shape", "Define a shape by the decimal of a bitvector, where 0 symbolizes a "
                      "position to be ignored, 1 a position considered. Default: ungapped.");
    parser.add_option(args.seed, 's', "seed", "Define seed.");

    try
    {
        parser.parse();
    }
    catch (seqan3::parser_invalid_argument const & ext)
    {
        seqan3::debug_stream << "Error. Incorrect command line input for IBF construct." << ext.what() << "\n";
        return -1;
    }

    if (args.samples.empty()) // If no samples are given, every file is seen as on experiment
    {
        args.samples.assign(args.sequence_files.size(),1);
    }
    // If sum of args.samples is not equal to number of files
    else if (std::accumulate(args.samples.rbegin(), args.samples.rend(), 0) != args.sequence_files.size())
    {
        seqan3::debug_stream << "Error. Incorrect command line input for multiple-samples." << "\n";
        return -1;
    }

    ibf(args);

    return 0;
}*/
// =====================================================================================================================
// search
// =====================================================================================================================
int main(int argc, char const ** argv)//int run_needle_search(argument_parser & parser)
{
    seqan3::argument_parser parser("needle-search", argc, argv);
    cmd_arguments args{};
    parser.info.author = "Mitra Darvish";
    parser.info.short_description = "Search through an IBF.";
    parser.add_positional_option(args.gene_file, "Please provide a sequence file.");
    parser.add_option(args.exp_file, 'x', "expression_file", "A tab seperated file containing expression information "
                      "per transcript given. By Default: Search for Existence in experiments",
                      seqan3::option_spec::DEFAULT, seqan3::input_file_validator{{"tsv"}});
    parser.add_option(args.path_in, 'i', "in", "Directory where input files can be found.");
    parser.add_option(args.expression, 'e', "expression", "Which expression level should be considered during a "
                      "search.");
    parser.add_option(args.k, 'k', "kmer", "Define kmer size.");
    parser.add_option(args.window_size, 'w', "window", "Define window size.");
    parser.add_option(args.shape, 'p', "shape", "Define a shape by the decimal of a bitvector, where 0 symbolizes a "
                      "position to be ignored, 1 a position considered. Default: ungapped.");
    parser.add_option(args.seed, 's', "seed", "Define seed.");
    parser.add_flag(args.compressed, 'c', "compressed", "If set ibf is compressed. Default: Not compressed.");

    try
    {
        parser.parse();
    }
    catch (seqan3::parser_invalid_argument const & ext)                     // catch user errors
    {
        seqan3::debug_stream << "Error. Incorrect command line input for IBF construct." << ext.what() << "\n";
        return -1;
    }

    if (args.compressed)
    {
        seqan3::binning_directory_compressed bd;
        return search(bd, args);
    }
    else
    {
        seqan3::binning_directory bd;
        return search(bd, args);
    }

    return 0;
}
// =====================================================================================================================
// main
// =====================================================================================================================
/*int main(int argc, char const ** argv)
{
    seqan3::argument_parser ibf_parser("needle-IBF", argc, argv);
    seqan3::argument_parser search_parser("needle-search", argc, argv);
    cmd_arguments args{};

    argument_parser top_level_parser{"needle", argc, argv, true};
    // Add information and flags to your top-level parser just as you would do with a normal one.
    // Note that all flags directed at the top-level parser must be specified BEFORE the subcommand key word.
    // Because of ambiguity, we do not allow any (positional) options for the top-level parser.
    top_level_parser.info.description.push_back("You can create an ibf and search in an ibf.");
    try
    {
        top_level_parser.parse(); // trigger command line parsing
    }
    catch (parser_invalid_argument const & ext) // catch user errors
    {
        debug_stream << "Error. Incorrect command. Use either needle ibf or needle search." << ext.what() << "\n";
        return -1;
    }
    argument_parser & sub_parser = top_level_parser.get_sub_parser(); // hold a reference to the sub_parser
    std::cout << "Proceed to sub parser." << std::endl;
    if (sub_parser.info.app_name == std::string_view{"needle-ibf"})
        run_needle_ibf(sub_parser);
    else if (sub_parser.info.app_name == std::string_view{"needle-search"})
        run_needle_search(sub_parser);
    else
        throw std::logic_error{"I do not know sub parser " + sub_parser.info.app_name};
}*/
