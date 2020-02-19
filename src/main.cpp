#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>

#include "give_me_five.hpp"

int main(int argc, char ** argv)
{
    seqan3::argument_parser parser{"I-want-five", argc, argv};

    // Declarations for input variables
    std::vector<size_t> in;
    std::filesystem::path out;
    bool yes = false;

    // Parser
    parser.info.author = "SeqAn-Team"; // give parser some infos
    parser.info.version = "1.0.0";
    parser.add_option(in, 'i', "input", "Input Data");
    parser.add_option(out, 'o', "output", "Output path");
    parser.add_flag(yes, 'y', "yes", "Set yes to true"); // example for a flag

    try
    {
         parser.parse();                                                  // trigger command line parsing
    }
    catch (seqan3::argument_parser_error const & ext)                     // catch user errors
    {
        seqan3::debug_stream << "Parsing error. No fives can be given. " << ext.what() << "\n"; // give error message
        return -1;
    }

    seqan3::debug_stream << "Hello world\n";
    if (yes) // if flag is set
        seqan3::debug_stream << "My app returns " << my_app::give_me_five() << ".\n";

    return 0;
}

