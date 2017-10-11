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
#include "PbIndexDump.h"
#include "PbIndexDumpVersion.h"
#include "Settings.h"
#include <cassert>
#include <cstddef>
#include <iostream>
using namespace std;

static
pbindexdump::Settings fromCommandLine(optparse::OptionParser& parser,
                                      int argc,
                                      char* argv[])
{
    const optparse::Values options = parser.parse_args(argc, argv);
    pbindexdump::Settings settings;

    // input
    const vector<string> positionalArgs = parser.args();
    const size_t numPositionalArgs = positionalArgs.size();
    if (numPositionalArgs == 0)
        settings.inputPbiFilename_ = "-"; // stdin
    else if (numPositionalArgs == 1)
        settings.inputPbiFilename_ = parser.args().front();
    else {
        assert(numPositionalArgs > 1);
        settings.errors_.emplace_back("pbindexdump does not support more than one input file per run");
    }

    // output format
    if (options.is_set("format"))
        settings.format_ = options["format"];

    // JSON options
    if (settings.format_ == "json") {
        if (options.is_set("json_indent_level"))
            settings.jsonIndentLevel_ = options.get("json_indent_level");
        if (options.is_set("json_raw"))
            settings.jsonRaw_ = options.get("json_raw");
    } else {
        if (options.is_set("json_indent_level") ||
            options.is_set("json_raw"))
        {
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
    ioGroup.add_option("")
           .dest("input")
           .metavar("input")
           .help("Input PBI file. If not provided, stdin will be used as input.");
    ioGroup.add_option("--format")
           .dest("format")
           .metavar("STRING")
           .help("Output format, one of:\n"
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
             .help("Prints fields in a manner that more closely reflects the PBI"
                   " file format - presenting data as per-field columns, not"
                   " per-record objects.");
    parser.add_option_group(jsonGroup);

    // parse command line for settings
    const pbindexdump::Settings settings = fromCommandLine(parser, argc, argv);
    if (!settings.errors_.empty()) {
        cerr << endl;
        for (const auto e : settings.errors_)
            cerr << "ERROR: " << e << endl;
        cerr << endl;
        parser.print_help();
        return EXIT_FAILURE;
    }

    // run tool
    try {
        pbindexdump::PbIndexDump::Run(settings);
        return EXIT_SUCCESS;
    }
    catch (std::exception& e) {
        cerr << "ERROR: " << e.what() << endl;
        return EXIT_FAILURE;
    }
}
