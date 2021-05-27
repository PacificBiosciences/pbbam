#include "PbBamifySettings.h"

#include <stdexcept>

#include <pbcopper/data/CigarOperation.h>

#include "PbBamifyVersion.h"

namespace PacBio {
namespace PbBamify {
namespace Options {

// clang-format off
const CLI_v2::Option InputFile{
R"({
    "names" : ["input"],
    "description" : "The aligned non-PacBio BAM file. If not provided, stdin will be used as input.",
    "type" : "file",
    "default" : ""
})"};

const CLI_v2::Option OutputFile{
R"({
    "names" : ["output"],
    "description" : "Path to the output BAM file. If not specified, output will be to the stdout.",
    "type" : "file",
    "default" : ""
})"};

const CLI_v2::Option VerboseLevel{
R"({
    "names" : ["verbose-level"],
    "description" : [
        "Specifies the level of info which will be output produced on stderr. ",
        "0 turns all output off, 1 outputs only warnings, while levels 2 and ",
        "above outputs a status message every 1000000 (2), 100000 (3), 1000 (4), ",
        "100 (5), 10 (6) and 1 (7) reads."
    ],
    "type" : "int",
    "default" : 3
})"};

const CLI_v2::PositionalArgument ReferenceFile{
R"({
    "name" : "ref.fa",
    "description" : "Reference used to align the input.",
    "type" : "file"
})"};

const CLI_v2::PositionalArgument PbbamReadFile{
R"({
    "name" : "IN.bam",
    "description" : "Input file(s). Maybe one of: DataSetXML, BAM file(s), or FOFN",
    "type" : "file"
})"};

// clang-format on

}  // namespace Options

CLI_v2::Interface Settings::CreateCLI()
{
    // clang-format off
    const std::string description{
        "pbbamify converts an arbitray aligned BAM file to a PacBio-compatible BAM file."
        "Input BAM file is read from a file or stdin, the raw-reads PacBio BAM is given"
        "as a parameter, and BAM output is written to stdout."
    };

    CLI_v2::Interface interface{"pbbamify", description, PbBamify::Version};
    interface.DisableNumThreadsOption();

    interface.AddOptions({
        Options::InputFile,
        Options::OutputFile,
        Options::VerboseLevel
    });
    interface.AddPositionalArguments({
        Options::ReferenceFile,
        Options::PbbamReadFile
    });

    Logging::LogConfig logConfig{Logging::LogLevel::INFO};
    logConfig.Fields = Logging::LogField::TIMESTAMP | Logging::LogField::LOG_LEVEL;
    interface.LogConfig(logConfig);
    interface.DefaultLogLevel(Logging::LogLevel::INFO); // fixme
    // clang-format on

    return interface;
}

Settings::Settings(const CLI_v2::Results& args)
    : InputFilename(args[Options::InputFile])
    , OutputFilename(args[Options::OutputFile])
    , VerboseLevel{args[Options::VerboseLevel]}
{
    // Reference & unaligned PacBio BAM files
    const auto& posArgs = args.PositionalArguments();
    if (posArgs.size() != 2) {
        throw std::runtime_error{"exactly two positional arguments must be provided"};
    }
    ReferenceFilename = posArgs[0];
    PbbamFilename = posArgs[1];

    // Input non-PacBio BAM
    if (InputFilename.empty()) {
        InputFilename = "-";
    }

    // Output aligned PacBio BAM
    if (OutputFilename.empty()) {
        OutputFilename = "-";
    }

    // Verbosity
    if (VerboseLevel < 0) {
        VerboseLevel = 0;
    }

    // Allow 'M' tags
    Data::CigarOperation::DisableAutoValidation();
}

}  // namespace PbBamify
}  // namespace PacBio
