// Copyright (c) 2015, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

// Author: Derek Barnett

#include "../common/OptionParser.h"
#include "Bam2Sam.h"
#include "Bam2SamVersion.h"
#include <string>
#include <vector>
#include <cassert>
#include <cstddef>
#include <cstdlib>

static
bam2sam::Settings fromCommandLine(optparse::OptionParser& parser,
                                    int argc, char* argv[])
{
    bam2sam::Settings settings;

    const optparse::Values options = parser.parse_args(argc, argv);

    // input
    const std::vector<std::string> positionalArgs = parser.args();
    const size_t numPositionalArgs = positionalArgs.size();
    if (numPositionalArgs == 0)
        settings.inputFilename_ = "-"; // stdin
    else if (numPositionalArgs == 1)
        settings.inputFilename_ = parser.args().front();
    else {
        assert(numPositionalArgs > 1);
        settings.errors_.emplace_back("bam2sam does not support more than one input file per run");
    }

    // header options
    if (options.is_set("no_header"))
        settings.noHeader_ = options.get("no_header");
    if (options.is_set("header_only"))
        settings.printHeaderOnly_ = options.get("header_only");

    if (settings.noHeader_ && settings.printHeaderOnly_)
        settings.errors_.emplace_back("conflicting arguments requested: --no-header and --header-only");

    return settings;
}

int main(int argc, char* argv[])
{
    // setup help & options
    optparse::OptionParser parser;
    parser.description("bam2sam converts a BAM file to SAM. It is essentially a stripped-down "
                       "'samtools view', mostly useful for testing/debugging without requiring samtools. "
                       "Input BAM file is read from a file or stdin, and SAM output is written to stdout."
                       );
    parser.prog("bam2sam");
    parser.usage("bam2sam [options] [input]");
    parser.version(bam2sam::Version);
    parser.add_version_option(true);
    parser.add_help_option(true);

    auto optionGroup = optparse::OptionGroup(parser, "Options");
    optionGroup.add_option("")
               .dest("input")
               .metavar("input")
               .help("Input BAM file. If not provided, stdin will be used as input.");
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
    }
    catch (std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
}
