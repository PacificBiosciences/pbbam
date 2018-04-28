// Author: Ivan Sovic (based on code from Derek Barnett)

#include <cassert>
#include <iostream>
#include <istream>

#include <pbbam/BamReader.h>
#include <pbbam/BamWriter.h>
#include <pbbam/DataSet.h>
#include <pbbam/FastaReader.h>
#include "../common/OptionParser.h"
#include "PbBamify.h"
#include "PbBamifyVersion.h"
#include "QueryLookup.h"

namespace PacBio {
namespace BAM {
namespace pbbamify {

class Settings
{
public:
    static Settings FromCommandLine(optparse::OptionParser& parser, int argc, char* argv[])
    {
        pbbamify::Settings settings;
        const optparse::Values options = parser.parse_args(argc, argv);

        // Pbbam file for BAM tags.
        const std::vector<std::string> positionalArgs = parser.args();
        if (positionalArgs.size() != 2) {
            settings.errors_.push_back("Exactly two positional arguments must be specified.");
        } else {
            settings.referenceFilename_ = positionalArgs[0];
            settings.pbbamFilename_ = positionalArgs[1];
        }

        // Input generic BAM file to turn into Pbbam. Optional, so
        // that BAM records can be piped in and conversion made
        // on the fly.
        if (options.is_set("input")) {
            settings.inputFilename_ = options["input"];
        } else {
            settings.inputFilename_ = "-";
        }

        // If not specified, output is to stdout.
        if (options.is_set("output")) {
            settings.outputFilename_ = options["output"];
        } else {
            settings.outputFilename_ = "-";
        }

        // Info messages can be written to stderr if verbose_level > 0.
        std::istringstream iss(std::string{options["verbose_level"]});
        iss >> settings.verboseLevel_;
        if (settings.verboseLevel_ < 0) {
            settings.verboseLevel_ = 0;
        }

        // Disable validation of CIGARs that might contain 'M'
        CigarOperation::validate_ = false;

        return settings;
    }

public:
    std::string inputFilename_;
    std::string outputFilename_;
    std::string referenceFilename_;
    std::string pbbamFilename_;
    std::vector<std::string> errors_;
    int32_t verboseLevel_;

private:
    Settings() {}
};

}  // namespace pbbamify
}  // namespace BAM
}  // namespace PacBio

int main(int argc, char* argv[])
{
    // Setup help & options
    optparse::OptionParser parser;
    parser.description(
        "pbbamify converts an arbitray aligned BAM file to a PacBio-compatible BAM file."
        "Input BAM file is read from a file or stdin, the raw-reads PacBio BAM is given"
        "as a parameter, and BAM output is written to stdout.");
    parser.prog("pbbamify");
    parser.usage("pbbamify [options] <ref.fa> <pb.bam>|<pb.fofn>|<pb.xml>");
    parser.version(PacBio::BAM::pbbamify::Version);
    parser.add_version_option(true);
    parser.add_help_option(true);

    parser.set_defaults("verbose_level", "3");

    auto optionGroup = optparse::OptionGroup(parser, "Options");
    optionGroup.add_option("").dest("ref").help("Reference used to align the input.");
    optionGroup.add_option("--input").dest("input").metavar("STR").help(
        "The aligned non-PacBio BAM file. If not provided, stdin will be used as input.");
    optionGroup.add_option("--output")
        .dest("output")
        .metavar("STR")
        .help("Path to the output BAM file. If not specified, output will be to the stdout.");
    optionGroup.add_option("--verbose-level")
        .dest("verbose_level")
        .type("int")
        .metavar("INT")
        .set_default("3")
        .help(
            "Specifies the level of info which will be output produced on"
            "stderr. 0 turns all output off, 1 outputs only warnings, "
            "while levels 2 and above outputs a status message every "
            "1000000 (2), 100000 (3), 1000 (4), 100 (5), 10 (6) and 1 (7) reads.");
    optionGroup.add_option("").dest("pbbam").help("A PacBio BAM file containing raw reads.");
    // A Pbbam can be one of the following:
    // - DataSetXML
    // - FOFN
    // - BAM
    parser.add_option_group(optionGroup);

    // Parse command line for settingas.
    const PacBio::BAM::pbbamify::Settings settings =
        PacBio::BAM::pbbamify::Settings::FromCommandLine(parser, argc, argv);
    if (!settings.errors_.empty()) {
        std::cerr << std::endl;
        for (const auto& e : settings.errors_) {
            std::cerr << "ERROR: " << e << std::endl;
        }
        std::cerr << std::endl;
        parser.print_help();
        return EXIT_FAILURE;
    }

    // Run the tool.
    try {
        // setup our @PG entry to add to header
        PacBio::BAM::ProgramInfo pbbamifyProgram;
        pbbamifyProgram.Id(std::string("pbbamify-") + PacBio::BAM::pbbamify::Version)
            .Name("pbbamify")
            .Version(PacBio::BAM::pbbamify::Version);

        PacBio::BAM::DataSet dataset = PacBio::BAM::DataSet(settings.pbbamFilename_);
        PacBio::BAM::BamReader inputBamReader(settings.inputFilename_);
        PacBio::BAM::BamHeader newHeader;

        {  // A separate block to close the reference file after the header is formed.
            // Using a sequential reader to construct the header SN lines in order, fast.
            PacBio::BAM::FastaReader ref_reader(settings.referenceFilename_);
            newHeader =
                PacBio::BAM::pbbamify::Pbbamify::ComposeHeader(dataset, ref_reader, inputBamReader);
        }

        std::shared_ptr<PacBio::BAM::pbbamify::QueryLookup> queryLookup =
            PacBio::BAM::pbbamify::CreateQueryLookup(dataset);
        queryLookup->Load();

        {  // A block is used here to close the bamWriter and the reference reader.
            // (Even though this will be done as soon as the 'try' block ends, this safeguards if any
            // code should be added in between at some point.)
            PacBio::BAM::IndexedFastaReader indexedRefReader(settings.referenceFilename_);
            PacBio::BAM::BamWriter bamWriter(settings.outputFilename_, newHeader);
            bool augment_rv = PacBio::BAM::pbbamify::Pbbamify::AugmentAlignments(
                queryLookup, indexedRefReader, inputBamReader, bamWriter, settings.verboseLevel_);
            if (augment_rv == false) {
                return EXIT_FAILURE;
            }
        }

        return EXIT_SUCCESS;
    } catch (std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
}
