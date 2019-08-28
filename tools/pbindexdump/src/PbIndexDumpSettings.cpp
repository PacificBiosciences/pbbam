// Author: Derek Barnett

#include "PbIndexDumpSettings.h"

#include <stdexcept>

#include "PbIndexDumpVersion.h"

namespace PacBio {
namespace PbIndexDump {
namespace Options {

// clang-format off
const CLI_v2::PositionalArgument InputFile{
R"({
    "name" : "input.bam.pbi",
    "description" : "Input PBI file. If not provided, stdin will be used as input.",
    "type" : "file",
    "required" : false
})"};

const CLI_v2::Option Format{
R"({
    "names" : ["format"],
    "description" : "Output format.",
    "type" : "string",
    "choices" : ["json", "cpp"],
    "default" : "json"
})"};

const CLI_v2::Option JsonIndentLevel{
R"({
    "names" : ["json-indent-level"],
    "description" : "JSON indent level.",
    "type" : "int",
    "default" : 4
})"};

const CLI_v2::Option JsonRaw{
R"({
    "names" : ["json-raw"],
    "description" : [
        "Print fields in a layout that more closely reflects the PBI ",
        "file format (per-field columns, not per-record objects."
    ]
})"};

// clang-format on

}  // namespace Options

CLI_v2::Interface Settings::CreateCLI()
{
    // clang-format off
    const std::string description{
        "pbindexdump prints a human-readable view of PBI data to stdout."
    };

    CLI_v2::Interface interface{"pbindexdump", description, PbIndexDump::Version};
    interface.DisableLogFileOption()
             .DisableLogLevelOption()
             .DisableNumThreadsOption();

    interface.AddPositionalArguments({
        Options::InputFile
    });

    interface.AddOptionGroup("Output Options", {
        Options::Format,
        Options::JsonIndentLevel,
        Options::JsonRaw
    });

    interface.HelpFooter({
        "Supported output formats:\n"
        "  json: 'pretty-printed' JSON\n"
        "  cpp:  copy/paste-able C++ code that can be used to construct the\n"
        "        equivalent BAM::PbiRawData object."
    });
    // clang-format on

    return interface;
}

Settings::Settings(const CLI_v2::Results& args)
    : Format(args[Options::Format])
    , JsonIndentLevel{args[Options::JsonIndentLevel]}
    , JsonRaw(args[Options::JsonRaw])
{
    // input file
    const auto& posArgs = args.PositionalArguments();
    if (posArgs.empty())
        InputFile = "-";
    else if (posArgs.size() == 1)
        InputFile = posArgs[0];
    else
        throw std::runtime_error{"too many arguments provided."};

    // format sanity check
    if (Format != "json" && Format != "cpp")
        throw std::runtime_error{"unsupported format requested: '" + Format + "'"};

    // JSON options sanity check
    if (Format != "json") {
        if (args[Options::JsonRaw].IsUserProvided() ||
            args[Options::JsonIndentLevel].IsUserProvided()) {
            throw std::runtime_error{"JSON formatting options are not valid on non-JSON output"};
        }
    }
}

}  // namespace PbIndexDump
}  // namespace PacBio
