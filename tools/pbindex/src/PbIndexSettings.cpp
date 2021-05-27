#include "PbIndexSettings.h"

#include <stdexcept>

#include "PbIndexVersion.h"

namespace PacBio {
namespace PbIndex {
namespace Options {

// clang-format off
const CLI_v2::PositionalArgument InputFile{
R"({
    "name" : "IN.bam",
    "description" : "Input BAM file",
    "type" : "file"
})"};
// clang-format on

}  // namespace Options

PacBio::CLI_v2::Interface Settings::CreateCLI()
{
    // clang-format off
    const std::string description{
        "pbindex creates a index file that enables random-access to PacBio-specific "
        "data in BAM files. Generated index filename will be the same as input BAM "
        "plus .pbi suffix."
    };

    CLI_v2::Interface interface{"pbindex", description, PbIndex::Version};
    interface.DisableLogFileOption()
             .DisableLogLevelOption()
             .DisableNumThreadsOption();

    interface.AddPositionalArguments({
        Options::InputFile
    });
    // clang-format on

    return interface;
}

Settings::Settings(const CLI_v2::Results& args) : InputFile(args[Options::InputFile]) {}

}  // namespace PbIndex
}  // namespace PacBio
