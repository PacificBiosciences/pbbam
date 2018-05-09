// Author: Derek Barnett

#include <cassert>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include "../common/OptionParser.h"

#include "PbIndexDump.h"
#include "PbIndexDumpVersion.h"
#include "Settings.h"

static pbindexdump::Settings fromCommandLine(optparse::OptionParser& parser, int argc, char* argv[])
{
    const optparse::Values options = parser.parse_args(argc, argv);
    pbindexdump::Settings settings;

    // input
    const std::vector<std::string> positionalArgs = parser.args();
    const size_t numPositionalArgs = positionalArgs.size();
    if (numPositionalArgs == 0)
        settings.inputPbiFilename_ = "-";  // stdin
    else if (numPositionalArgs == 1)
        settings.inputPbiFilename_ = parser.args().front();
    else {
        assert(numPositionalArgs > 1);
        settings.errors_.emplace_back(
            "pbindexdump does not support more than one input file per run");
    }

    // output format
    if (options.is_set("format")) settings.format_ = options["format"];

    // JSON options
    if (settings.format_ == "json") {
        if (options.is_set("json_indent_level"))
            settings.jsonIndentLevel_ = options.get("json_indent_level");
        if (options.is_set("json_raw")) settings.jsonRaw_ = options.get("json_raw");
    } else {
        if (options.is_set("json_indent_level") || options.is_set("json_raw")) {
            settings.errors_.emplace_back("JSON formatting options not valid on non-JSON output");
        }
    }

    return settings;
}

int main(int argc, char* argv[])
{
    // setup help & options
    optparse::OptionParser parser;
    parser.description("pbindexdump prints a human-readable view of PBI data to stdout.");
    parser.prog("pbindexdump");
    parser.usage("pbindexdump [options] [input]");
    parser.version(pbindexdump::Version);
    parser.add_version_option(true);
    parser.add_help_option(true);

    auto ioGroup = optparse::OptionGroup(parser, "Input/Output");
    ioGroup.add_option("").dest("input").metavar("input").help(
        "Input PBI file. If not provided, stdin will be used as input.");
    ioGroup.add_option("--format")
        .dest("format")
        .metavar("STRING")
        .help(
            "Output format, one of:\n"
            "    json, cpp\n\n"
            "json: pretty-printed JSON [default]\n\n"
            "cpp: copy/paste-able C++ code that can be used to construct the"
            " equivalent PacBio::BAM::PbiRawData object");
    parser.add_option_group(ioGroup);

    auto jsonGroup = optparse::OptionGroup(parser, "JSON Formatting");
    jsonGroup.add_option("--json-indent-level")
        .dest("json_indent_level")
        .metavar("INT")
        .help("JSON indent level [4]");
    jsonGroup.add_option("--json-raw")
        .dest("json_raw")
        .action("store_true")
        .help(
            "Prints fields in a manner that more closely reflects the PBI"
            " file format - presenting data as per-field columns, not"
            " per-record objects.");
    parser.add_option_group(jsonGroup);

    // parse command line for settings
    const pbindexdump::Settings settings = fromCommandLine(parser, argc, argv);
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
        pbindexdump::PbIndexDump::Run(settings);
        return EXIT_SUCCESS;
    } catch (std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
}
