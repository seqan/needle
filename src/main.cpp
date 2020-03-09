#include <sstream>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>

#include "fastq_conversion.hpp"

int main(int argc, char ** argv)
{
    seqan3::argument_parser parser{"Fastq-to-Fasta-Converter", argc, argv};

    // Declarations for argument parser
    std::filesystem::path fastq_file{};
    std::filesystem::path output_file{};
    bool verbose = false;

    // Parser
    parser.info.author = "SeqAn-Team"; // give parser some infos
    parser.info.version = "1.0.0";
    parser.add_positional_option(fastq_file, "Please provide a fastq file.",
                                 seqan3::input_file_validator{{"fq","fastq"}}); // Takes a fastq file and validates it
    //output path as option, otherwise output is printed
    parser.add_option(output_file, 'o', "output", "The file for fasta output. Default: stdout");
    parser.add_flag(verbose, 'v', "verbose", "Give more detailed information here."); // example for a flag

    try
    {
         parser.parse();                                                  // trigger command line parsing
    }
    catch (seqan3::argument_parser_error const & ext)                     // catch user errors
    {
        seqan3::debug_stream << "Parsing error. " << ext.what() << "\n"; // give error message
        return -1;
    }

    convert_fastq(fastq_file, output_file); // Call fastq to fasta converter

    if (verbose) // if flag is set
        seqan3::debug_stream << "Conversion was a sucess. Congrats!\n";



    return 0;
}
