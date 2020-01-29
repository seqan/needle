#include <seqan3/core/debug_stream.hpp>

#include "give_me_five.hpp"

int main()
{
    seqan3::debug_stream << "Hello world\n";
    seqan3::debug_stream << "My app returns " << my_app::give_me_five() << ".\n";
    return 0;
}

