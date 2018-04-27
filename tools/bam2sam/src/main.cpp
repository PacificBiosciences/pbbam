// Author: Derek Barnett

#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <string>
#include <vector>
#include "../common/OptionParser.h"
#include "Bam2Sam.h"
#include "Bam2SamVersion.h"

static bam2sam::Settings fromCommandLine(optparse::OptionParser& parser, int argc, char* argv[])
{
    bam2sam::Settings settings;

    const optparse::Values options = parser.parse_args(argc, argv);

    // input
    const std::vector<std::string> positionalArgs = parser.args();
    const size_t numPositionalArgs = positionalArgs.size();
    if (numPositionalArgs == 0)
        settings.inputFilename_ = "-";  // stdin
    else if (numPositionalArgs == 1)
        settings.inputFilename_ = parser.args().front();
    else {
        assert(numPositionalArgs > 1);
        settings.errors_.emplace_back("bam2sam does not support more than one input file per run");
    }

    // header options
    if (options.is_set("no_header")) settings.noHeader_ = options.get("no_header");
    if (options.is_set("header_only")) settings.printHeaderOnly_ = options.get("header_only");

    if (settings.noHeader_ && settings.printHeaderOnly_)
        settings.errors_.emplace_back(
            "conflicting arguments requested: --no-header and --header-only");

    return settings;
}

int main(int argc, char* argv[])
{
    // setup help & options
    optparse::OptionParser parser;
    parser.description(
        "bam2sam converts a BAM file to SAM. It is essentially a stripped-down "
        "'samtools view', mostly useful for testing/debugging without requiring samtools. "
        "Input BAM file is read from a file or stdin, and SAM output is written to stdout.");
    parser.prog("bam2sam");
    parser.usage("bam2sam [options] [input]");
    parser.version(bam2sam::Version);
    parser.add_version_option(true);
    parser.add_help_option(true);

    auto optionGroup = optparse::OptionGroup(parser, "Options");
    optionGroup.add_option("").dest("input").metavar("input").help(
        "Input BAM file. If not provided, stdin will be used as input.");
    optionGroup.add_option("--no-header")
        .dest("no_header")
        .action("store_true")
        .help("Omit header from output.");
    optionGroup.add_option("--header-only")
        .dest("header_only")
        .action("store_true")
        .help("Print only the header (no records).");
    parser.add_option_group(optionGroup);

    // parse command line for settings
    const bam2sam::Settings settings = fromCommandLine(parser, argc, argv);
    if (!settings.errors_.empty()) {
        std::cerr << std::endl;
        for (const auto e : settings.errors_)
            std::cerr << "ERROR: " << e << std::endl;
        std::cerr << std::endl;
        parser.print_help();
        return EXIT_FAILURE;
    }

    // run tool
    try {
        bam2sam::PbBam2Sam::Run(settings);
        return EXIT_SUCCESS;
    } catch (std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
}
