// Author: Derek Barnett

#include <cassert>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include "../common/OptionParser.h"
#include "PbIndex.h"
#include "PbIndexVersion.h"

static pbindex::Settings fromCommandLine(optparse::OptionParser& parser, int argc, char* argv[])
{
    const optparse::Values options = parser.parse_args(argc, argv);
    //    ()options;

    pbindex::Settings settings;

    // get input filename
    const std::vector<std::string> positionalArgs = parser.args();
    const size_t numPositionalArgs = positionalArgs.size();
    if (numPositionalArgs == 0)
        settings.errors_.push_back("pbindex requires an input BAM filename");
    else if (numPositionalArgs == 1)
        settings.inputBamFilename_ = parser.args().front();
    else {
        assert(numPositionalArgs > 1);
        settings.errors_.push_back("pbindex does not support more than one input file per run");
    }

    return settings;
}

int main(int argc, char* argv[])
{
    // setup help & options
    optparse::OptionParser parser;
    parser.description(
        "pbindex creates a index file that enables random-access to PacBio-specific data in BAM "
        "files. "
        "Generated index filename will be the same as input BAM plus .pbi suffix.");
    parser.prog("pbindex");
    parser.usage("pbindex <input>");
    parser.version(pbindex::Version);
    parser.add_version_option(true);
    parser.add_help_option(true);

    auto ioGroup = optparse::OptionGroup(parser, "Input/Output");
    ioGroup.add_option("").dest("input").metavar("input").help("Input BAM file");
    parser.add_option_group(ioGroup);

    // parse command line for settings
    const pbindex::Settings settings = fromCommandLine(parser, argc, argv);
    if (!settings.errors_.empty()) {
        std::cerr << std::endl;
        for (const auto e : settings.errors_)
            std::cerr << "ERROR: " << e << std::endl;
        std::cerr << std::endl;
        parser.print_help();
        return EXIT_FAILURE;
    }

    // run tool
    return pbindex::PbIndex::Run(settings);
}
