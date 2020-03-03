#include <seqan3/argument_parser/all.hpp>
#include "ibf.h"
#include "search.h"

void initialise_argument_parser(seqan3::argument_parser & parser, arguments & args)
{
    parser.add_flag(args.compressed, 'c', "compressed", "If c is set, ibf is compressed. Default: Not compressed.");
    parser.add_option(args.k, 'k', "kmer", "Define kmer size.");
    parser.add_option(args.window_size, 'w', "window", "Define window size.");
    parser.add_option(args.shape, 'p', "shape", "Define a shape by the decimal of a bitvector, where 0 symbolizes a "
                                                "position to be ignored, 1 a position considered. Default: ungapped.");
    parser.add_option(args.seed, 's', "seed", "Define seed.");
}

int run_needle_ibf(seqan3::argument_parser & parser)
{
    arguments args{};
    ibf_arguments ibf_args{};
    parser.info.short_description = "Constructs an IBF.";
    parser.info.version = "1.0.0";
    parser.info.author = "Mitra Darvish";

    parser.add_positional_option(ibf_args.sequence_files, "Please provide at least one sequence file.");
    parser.add_option(ibf_args.genome_file, 'g', "genom-mask", "Genom file used as a mask.");
    parser.add_option(ibf_args.bin_size, 'b', "bin-size", "List of bin sizes per expression level. If only one is given"
                                                          ", then that bin size is used for all expression levels.");
    parser.add_option(ibf_args.num_hash, 'n', "hash", "Number of hash functions that should be used when constructing "
                      "one IBF.");
    parser.add_option(ibf_args.path_out, 'o', "out", "Directory, where output files should be saved.");
    parser.add_option(ibf_args.expression_levels, 'e', "expression_levels", "Which expression levels should be used for"
                                                                            " constructing the IBFs.");
    parser.add_option(ibf_args.samples, 'm', "multiple-samples", "Define which samples belong together, sum has to be "
                                                                 "equal to number of sequence files. Default: Every"
                                                                 " sequence file is one sample from one experiment.");
    parser.add_option(ibf_args.cutoffs, 'u', "cut-offs", "Define for each sample, what number of found minimizers "
                                                         "should be considered the result of a sequencing error and "
                                                         "therefore be ignored. Default: Every sample has a cut off of "
                                                         "zero.");
    parser.add_option(ibf_args.normalization_method, 'a', "normalization-method", "Choose a normalization method: mean,"
                                                                                  " median or random. Default: median.");
    parser.add_option(ibf_args.random, 'r', "random-samples", "Choose the number of random sequences to pick from when "
                                                              "using normalization method random. Default: 1000.");
    parser.add_option(ibf_args.random, 'f', "experiment-names", "If set, names of the experiments are stored in a txt file.");

    initialise_argument_parser(parser, args);

    try
    {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        seqan3::debug_stream << "Error. Incorrect command line input for IBF construct." << ext.what() << "\n";
        return -1;
    }
    try
    {
        ibf(args, ibf_args);
    }
    catch (const std::invalid_argument & e)
    {
        std::cerr << e.what() << std::endl;
        return -1;
    }

    return 0;
}

int run_needle_search(seqan3::argument_parser & parser)
{
    arguments args{};
    search_arguments search_args{};
    parser.info.short_description = "Search through an IBF.";
    parser.info.version = "1.0.0";
    parser.info.author = "Mitra Darvish";

    parser.add_positional_option(search_args.search_file, "Please provide a sequence file.");
    parser.add_option(search_args.exp_file, 'x', "expression_file", "A tab seperated file containing expression "
                                                                    "information per transcript given. By Default: "
                                                                    "Search for Existence in experiments",
                      seqan3::option_spec::DEFAULT, seqan3::input_file_validator{{"tsv"}});
    parser.add_option(search_args.path_in, 'i', "in", "Directory where input files can be found.");
    parser.add_option(search_args.expression, 'e', "expression", "Which expression level should be considered during a "
                                                                 "search.");

    initialise_argument_parser(parser, args);

    try
    {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const & ext)                     // catch user errors
    {
        seqan3::debug_stream << "Error. Incorrect command line input for IBF construct." << ext.what() << "\n";
        return -1;
    }

    std::vector<uint32_t> results;
    try
    {
        results = search(args, search_args);
    }
    catch (const std::invalid_argument & e)
    {
        std::cerr << e.what() << std::endl;
        return -1;
    }

    std::cout << "Results:\n";
    for( auto & elem : results)
        std::cout << elem << " ";
    std::cout << "\n";
    return 0;
}

int main(int argc, char const ** argv)
{
    seqan3::argument_parser needle_parser{"needle", argc, argv, true, {"ibf", "search"}};
    needle_parser.info.description.push_back("Needle allows you to build an Interleaved Bloom Filter (IBF) with the "
                                             "command ibf or search an IBF with the search command.");
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
    if (sub_parser.info.app_name == std::string_view{"needle-ibf"})
        run_needle_ibf(sub_parser);
    else if (sub_parser.info.app_name == std::string_view{"needle-search"})
        run_needle_search(sub_parser);
    else
        throw std::logic_error{"The used sub parser is not known: " + sub_parser.info.app_name};
}
