#include <seqan3/argument_parser/all.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/core/debug_stream.hpp>

#include "minimiser.h"
#include "ibf.h"
#include "search.h"

uint32_t w_size;
uint64_t shape{};
uint64_t se;

void initialise_argument_parser(seqan3::argument_parser & parser, arguments & args)
{
    parser.add_flag(args.compressed, 'c', "compressed", "If c is set, ibf is compressed. Default: Not compressed.");
    parser.add_option(args.k, 'k', "kmer", "Define kmer size.");
    parser.add_option(w_size, 'w', "window", "Define window size. Default: 60.");
    parser.add_option(shape, 'p', "shape", "Define a shape by the decimal of a bitvector, where 0 symbolizes a "
                                           "position to be ignored, 1 a position considered. Default: ungapped.");
    parser.add_option(se, 's', "seed", "Define seed.");
}

void parsing(seqan3::argument_parser & parser, arguments & args)
{
    w_size = args.w_size.get();
    se = args.s.get();
    parser.parse();
    args.w_size = seqan3::window_size{w_size};
    if (shape == 0)
            args.shape = seqan3::ungapped{args.k};
    else
            args.shape = seqan3::bin_literal{shape};
    args.s = seqan3::seed{se};
}

// Initialize arguments for ibf and minimiser
void initialise_ibf_argument_parser(seqan3::argument_parser & parser, ibf_arguments & ibf_args)
{
    parser.info.version = "1.0.0";
    parser.info.author = "Mitra Darvish";

    parser.add_positional_option(ibf_args.sequence_files, "Please provide at least one sequence file.");
    parser.add_option(ibf_args.include_file, 'g', "genom-mask", "Genom file used as a mask.");
    parser.add_option(ibf_args.path_out, 'o', "out", "Directory, where output files should be saved.");
    parser.add_option(ibf_args.samples, 'm', "multiple-samples", "Define which samples belong together, sum has to be "
                                                                 "equal to number of sequence files. Default: Every"
                                                                 " sequence file is one sample from one experiment.");
    parser.add_flag(ibf_args.paired, 'q', "paired", "If set, experiments are paired. Default: Not paired.");
    parser.add_option(ibf_args.cutoffs, 'u', "cut-offs", "Define for each sample, what number of found minimisers "
                                                         "should be considered the result of a sequencing error and "
                                                         "therefore be ignored. Default: Every sample has a cut off of "
                                                         "zero.");
    parser.add_option(ibf_args.normalization_method, 'a', "normalization-method", "Choose a normalization method: mean,"
                                                                                  " median or random. Default: median.");
}

int run_needle_estimate(seqan3::argument_parser & parser)
{
    arguments args{};
    search_arguments search_args{};
    parser.info.short_description = "Estimate expression value of transcript based on IBFs.";
    parser.info.version = "1.0.0";
    parser.info.author = "Mitra Darvish";
    std::filesystem::path search_file;
    std::filesystem::path path_in{"./"};
    std::filesystem::path file_out{"expressions.out"};
    std::vector<uint64_t> expressions{};

    parser.add_positional_option(search_file, "Please provide a sequence file.");
    parser.add_option(file_out, 'o', "out", "File where output should be stored.");
    parser.add_option(expressions, 'e', "expression", "Which expression levels should be considered during a "
                                                      "search.");
    parser.add_option(path_in, 'i', "in", "Directory where input files can be found.");
    parser.add_option(search_args.threshold, 't', "threshold", "The minimal amount of minimisers found in a transcript"
                                                                " to consider it as found in an IBF. Default: 0.5");
    initialise_argument_parser(parser, args);

    try
    {
        parsing(parser, args);
    }
    catch (seqan3::argument_parser_error const & ext)                     // catch user errors
    {
        seqan3::debug_stream << "Error. Incorrect command line input for search. " << ext.what() << "\n";
        return -1;
    }

    try
    {
        call_estimate(args, search_args, expressions, file_out, search_file, path_in);
    }
    catch (const std::invalid_argument & e)
    {
        std::cerr << e.what() << std::endl;
        return -1;
    }

    return 0;
}


int run_needle_ibf(seqan3::argument_parser & parser)
{
    arguments args;
    initialise_argument_parser(parser, args);
    ibf_arguments ibf_args{};
    initialise_ibf_argument_parser(parser, ibf_args);

    parser.info.short_description = "Constructs an IBF.";
    parser.add_option(ibf_args.bin_size, 'b', "bin-size", "List of bin sizes per expression level. If only one is given"
                                                          ", then that bin size is used for all expression levels.");
    parser.add_option(ibf_args.expression_levels, 'e', "expression_levels", "Which expression levels should be used for"
                                                                            " constructing the IBFs. Default: [0.5,1,2,4].");
    parser.add_option(ibf_args.num_hash, 'n', "hash", "Number of hash functions that should be used when constructing "
                                                      "one IBF.");
    parser.add_option(ibf_args.experiment_names, 'f', "experiment-names", "If set, names of the experiments are stored"
                                                                          " in a txt file.");

    try
    {
        parsing(parser, args);
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        seqan3::debug_stream << "Error. Incorrect command line input for ibf. " << ext.what() << "\n";
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

int run_needle_ibf_min(seqan3::argument_parser & parser)
{
    arguments args{};
    ibf_arguments ibf_args{};
    std::vector<std::filesystem::path> minimiser_files{};
    std::filesystem::path header_file = ""; // if only one header file should be used
    float fpr{0.05}; // False Positive Rate

    parser.info.short_description = "Constructs an IBF from the minimiser and header files created by needle minimiser.";
    parser.add_flag(args.compressed, 'c', "compressed", "If c is set, ibf is compressed. Default: Not compressed.");
    parser.add_positional_option(minimiser_files, "Please provide at least one minimiser file. It is assumed that the "
                                                  "header file exits in the same directory.");
    parser.add_option(fpr, 'f', "fpr", "False positive rate for the IBF. Default: 0.05.");
    parser.add_option(ibf_args.expression_levels, 'e', "expression_levels", "Which expression levels should be used for"
                                                                            " constructing the IBFs. Default: The "
                                                                            "expression levels found in the header files.");
    parser.add_option(ibf_args.include_file, 'g', "genom-mask", "Genom file used as a mask.");
    parser.add_option(ibf_args.path_out, 'o', "out", "Directory, where output files should be saved.");
    parser.add_option(ibf_args.num_hash, 'n', "hash", "Number of hash functions that should be used when constructing "
                                                      "one IBF.");
    parser.add_option(ibf_args.normalization_method, 'a', "normalization-method", "Choose a normalization method: mean,"
                                                                                  " median or random. Default: median.");

    try
    {
        parsing(parser, args);
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        seqan3::debug_stream << "Error. Incorrect command line input for ibfmin. " << ext.what() << "\n";
        return -1;
    }

    try
    {
        ibf(minimiser_files, header_file, args, ibf_args, fpr);
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
    arguments args{};
    initialise_argument_parser(parser, args);
    ibf_arguments ibf_args{};
    initialise_ibf_argument_parser(parser, ibf_args);
    parser.info.short_description = "Calculates minimiser for given experiments.";
    parser.add_option(ibf_args.expression_levels, 'e', "expression_levels", "The expression levels are used for counting"
                                                                            " how many minimisers are greater or equal "
                                                                            "than the given expression levels. Default "
                                                                            "is an expression level of zero, so all "
                                                                            "minimisers are counted. Multiple levels can"
                                                                            "be given, so multiple counts will be "
                                                                            "calculated. Default: [0].");

    try
    {
        parsing(parser, args);
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        seqan3::debug_stream << "Error. Incorrect command line input for minimiser. " << ext.what() << "\n";
        return -1;
    }
    try
    {
        minimiser(args, ibf_args);
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
    parser.add_option(search_args.threshold, 't', "threshold", "The minimal amount of minimisers found in a transcript"
                                                                " to consider it as found in an IBF. Default: 0.5");

    initialise_argument_parser(parser, args);

    try
    {
        parsing(parser, args);
    }
    catch (seqan3::argument_parser_error const & ext)                     // catch user errors
    {
        seqan3::debug_stream << "Error. Incorrect command line input for search. " << ext.what() << "\n";
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

int run_needle_stats(seqan3::argument_parser & parser)
{
    arguments args{};
    ibf_arguments ibf_args{};
    std::vector<std::filesystem::path> minimiser_files{};

    parser.info.short_description = "Get statistics from header files produced by needle minimiser.";
    parser.info.version = "1.0.0";
    parser.info.author = "Mitra Darvish";

    parser.add_positional_option(minimiser_files, "Please provide at least one header file.");

    try
    {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const & ext)                     // catch user errors
    {
        seqan3::debug_stream << "Error. Incorrect command line input for stats. " << ext.what() << "\n";
        return -1;
    }

    std::vector<std::tuple<std::vector<uint64_t>, std::vector<uint64_t>>> results;

    try
    {
        results = statistics(args, ibf_args, minimiser_files);
    }
    catch (const std::invalid_argument & e)
    {
        std::cerr << e.what() << std::endl;
        return -1;
    }

    for (unsigned i = 0; i < results.size(); ++i)
    {
        std::cout << "For expression level " << std::get<0>(results[i])[0] << ":\n";
        std::cout << "Average normalized expression value: " << std::get<0>(results[i])[1] <<"\n";
        std::cout << "Minimum of Counts: " << std::get<1>(results[i])[0] << "\n";
        std::cout << "Median of Counts: " << std::get<1>(results[i])[1] << "\n";
        std::cout << "Maximum of Counts: " << std::get<1>(results[i])[2] << "\n\n\n";
    }


    return 0;
}

int main(int argc, char const ** argv)
{
    seqan3::argument_parser needle_parser{"needle", argc, argv, true, {"estimate", "ibf", "ibfmin", "minimiser", "search",
                                                                       "stats"}};
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
    if (sub_parser.info.app_name == std::string_view{"needle-estimate"})
        run_needle_estimate(sub_parser);
    else if (sub_parser.info.app_name == std::string_view{"needle-ibf"})
        run_needle_ibf(sub_parser);
    else if (sub_parser.info.app_name == std::string_view{"needle-ibfmin"})
        run_needle_ibf_min(sub_parser);
    else if (sub_parser.info.app_name == std::string_view{"needle-minimiser"})
        run_needle_minimiser(sub_parser);
    else if (sub_parser.info.app_name == std::string_view{"needle-search"})
        run_needle_search(sub_parser);
    else if (sub_parser.info.app_name == std::string_view{"needle-stats"})
        run_needle_stats(sub_parser);
    else
        throw std::logic_error{"The used sub parser is not known: " + sub_parser.info.app_name};
}
