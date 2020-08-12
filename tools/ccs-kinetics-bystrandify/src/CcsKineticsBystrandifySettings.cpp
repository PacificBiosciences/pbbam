#include "CcsKineticsBystrandifySettings.h"

#include <algorithm>
#include <stdexcept>

#include "CcsKineticsBystrandifyVersion.h"

namespace PacBio {
namespace CcsKineticsBystrandify {
namespace Options {

// clang-format off
const CLI_v2::Option MinCoverage{
R"({
    "names" : ["min-coverage"],
    "description" : [
        "Specifies the minimum number of passes per strand (fn/rn) ",
        "for creating a strand-specific read."
    ],
    "type" : "int",
    "default" : 1
})"};

const CLI_v2::PositionalArgument InputFile{
R"({
    "name" : "IN.bam",
    "description" : "Input BAM file",
    "type" : "file"
})"};

const CLI_v2::PositionalArgument OutputFile{
R"({
    "name" : "OUT.bam",
    "description" : "Output BAM file",
    "type" : "file"
})"};
// clang-format on

}  // namespace Options

CLI_v2::Interface Settings::CreateCLI()
{
    const std::string description{
        "ccs-kinetics-bystrandify converts a BAM containing CCS-Kinetics tags to a pseudo-bystrand "
        "CCS BAM with pw/ip tags that can be used as a substitute for subreads in applications "
        "expecting such kinetic information."};

    CLI_v2::Interface interface{"ccs-kinetics-bystrandify", description,
                                CcsKineticsBystrandify::Version};
    interface.DisableNumThreadsOption();

    interface.AddOptions({
        Options::MinCoverage,
    });
    interface.AddPositionalArguments({
        Options::InputFile, Options::OutputFile,
    });

    Logging::LogConfig logConfig{Logging::LogLevel::INFO};
    logConfig.Fields = Logging::LogField::TIMESTAMP | Logging::LogField::LOG_LEVEL;
    interface.LogConfig(logConfig);
    interface.DefaultLogLevel(Logging::LogLevel::INFO);

    return interface;
}

Settings::Settings(const CLI_v2::Results& args)
    : MinCoverage{std::max<int32_t>(1, args[Options::MinCoverage])}
{
    // Reference & unaligned PacBio BAM files
    const auto& posArgs = args.PositionalArguments();
    if (posArgs.size() != 2)
        throw std::runtime_error{"exactly two positional arguments must be provided"};
    InputFilename = posArgs[0];
    OutputFilename = posArgs[1];
}

}  // namespace CcsKineticsBystrandify
}  // namespace PacBio
