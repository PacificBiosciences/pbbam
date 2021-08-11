#include "Bam2SamSettings.h"

#include <stdexcept>

#include "Bam2SamVersion.h"

namespace PacBio {
namespace Bam2Sam {
namespace Options {

// clang-format off
const CLI_v2::Option NoHeader{
R"({
    "names" : ["no-header"],
    "description" : "Omit header from output."
})"};

const CLI_v2::Option HeaderOnly{
R"({
    "names" : ["header-only"],
    "description" : "Print only the header (no records)."
})"};

const CLI_v2::PositionalArgument InputFile{
R"({
    "name" : "IN.bam",
    "description" : "Input BAM file. If not provided, stdin will be used as input.",
    "type" : "file",
    "required" : false
})"};
// clang-format on

}  // namespace Options

CLI_v2::Interface Settings::CreateCLI()
{
    // clang-format off
    const std::string description{
        "bam2sam converts a BAM file to SAM. It is essentially a stripped-down\n"
        "'samtools view', mostly useful for testing/debugging without requiring samtools.\n"
        "Input BAM file is read from a file or stdin, and SAM output is written to stdout."};

    CLI_v2::Interface interface{"bam2sam", description, Bam2Sam::Version};
    interface.DisableLogFileOption()
             .DisableLogLevelOption()
             .DisableNumThreadsOption();

    interface.AddOptionGroup("Options",
    {
        Options::NoHeader,
        Options::HeaderOnly
    });
    interface.AddPositionalArguments({
        Options::InputFile
    });
    // clang-format on

    return interface;
}

Settings::Settings(const CLI_v2::Results& args)
    : NoHeader{args[Options::NoHeader]}, HeaderOnly{args[Options::HeaderOnly]}
{
    // input file
    const auto& posArgs = args.PositionalArguments();
    if (posArgs.empty()) {
        InputFilename = "-";
    } else {
        InputFilename = posArgs.front();
    }

    // validate header print mode
    if (NoHeader && HeaderOnly) {
        throw std::runtime_error{
            "conflicting arguments requested '--no-header' and '--header-only'"};
    }
}

}  // namespace Bam2Sam
}  // namespace PacBio
