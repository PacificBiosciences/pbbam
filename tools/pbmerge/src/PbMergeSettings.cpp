// Author: Derek Barnett

#include "PbMergeSettings.h"

#include <stdexcept>

#include "PbMergeVersion.h"

namespace PacBio {
namespace PbMerge {
namespace Options {

// clang-format off
const CLI_v2::Option OutputFile{
R"({
    "names" : ["o"],
    "description" : "Output BAM filename. Writes to stdout if not provided.",
    "type" : "string",
    "default" : ""
})"};

const CLI_v2::Option NoPbi{
R"({
    "names" : ["no-pbi"],
    "description" : "Disables creation of PBI index file. PBI always disabled when writing to stdout."
})"};

const CLI_v2::PositionalArgument InputFiles{
R"({
    "name" : "INPUT",
    "description" : "Input file(s). Maybe one of: DataSetXML, BAM file(s), or FOFN",
    "type" : "file"
})"};
// clang-format on

}  // namespace Options

CLI_v2::Interface Settings::CreateCLI()
{
    // clang-format off
    const std::string description{
        "pbmerge merges PacBio BAM files. If the input is DataSetXML, any filters will be applied."
    };

    CLI_v2::Interface interface{"pbmerge", description, PbMerge::Version};
    interface.DisableLogFileOption()
             .DisableLogLevelOption()
             .DisableNumThreadsOption();

    interface.AddOptionGroup("Input/Output", {
        Options::OutputFile,
        Options::NoPbi
    });
    interface.AddPositionalArguments({
        Options::InputFiles
    });
    interface.HelpFooter(R"(Examples:
    $ pbmerge -o merged.bam data.subreadset.xml
    $ pbmerge -o merged.bam data_1.bam data_2.bam data_3.bam
    $ pbmerge -o merged.bam data_bams.fofn)");

    // clang-format on
    return interface;
}

Settings::Settings(const CLI_v2::Results& args) : OutputFile(args[Options::OutputFile])
{
    // input file(s)
    const auto& posArgs = args.PositionalArguments();
    if (posArgs.empty()) throw std::runtime_error{"at least input one file must be specified"};
    InputFiles = posArgs;

    // output (stdout?)
    if (OutputFile.empty()) OutputFile = "-";

    // create PBI?
    if (OutputFile == "-")
        CreatePbi = false;  // always skip PBI if writing to stdout
    else
        CreatePbi = !args[Options::NoPbi];  // create PBI unless requested
}

}  // namespace PbMerge
}  // namespace PacBio
