// Copyright (c) 2014, Pacific Biosciences of California, Inc.
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

#include "pbbam/SamHeaderCodec.h"
#include <vector>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

SamHeader SamHeaderCodec::Decode(const string& text)
{
    const string& Prefix_HD = string("@HD");
    const string& Prefix_SQ = string("@SQ");
    const string& Prefix_RG = string("@RG");
    const string& Prefix_PG = string("@PG");
    const string& Prefix_CO = string("@CO");

    SamHeader result;

    istringstream s(text);
    string line("");
    string firstToken;
    string restOfLine;
    while (getline(s, line)) {

        // skip if line is not long enough to contain true values
        if (line.length() < 5)
            continue;

        // determine token at beginning of line
        firstToken = line.substr(0,3);
        restOfLine = line.substr(4);
        if      (firstToken == Prefix_HD) DecodeHeaderLine(restOfLine, &result);
        else if (firstToken == Prefix_SQ) DecodeSequenceLine(restOfLine, &result);
        else if (firstToken == Prefix_RG) DecodeReadGroupLine(restOfLine, &result);
        else if (firstToken == Prefix_PG) DecodeProgramLine(restOfLine, &result);
        else if (firstToken == Prefix_CO) DecodeCommentLine(restOfLine, &result);
    }

    return result;
}

static
const vector<string> Split(const string& line, const char delim = '\t')
{
    vector<string> tokens;
    stringstream lineStream(line);
    string token;
    while (getline(lineStream, token, delim))
        tokens.push_back(token);
    return tokens;
}

void SamHeaderCodec::DecodeHeaderLine(const string& line,
                                      SamHeader* header)
{
    static string Token_VN = string("VN");
    static string Token_SO = string("SO");
    static string Token_pb = string("pb");

    // split HD lines into tokens
    const vector<string>& tokens = Split(line);

    // iterate over tokens
    auto tokenEnd  = tokens.cend();
    for (auto tokenIter = tokens.cbegin(); tokenIter != tokenEnd; ++tokenIter) {

        // get tag/value
        const string& tokenTag   = (*tokenIter).substr(0,2);
        const string& tokenValue = (*tokenIter).substr(3);

        // set header contents
        if      (tokenTag == Token_VN) header->version   = tokenValue;
        else if (tokenTag == Token_SO) header->sortOrder = tokenValue;
        else if (tokenTag == Token_pb) header->pacbioBamVersion = tokenValue;
    }

    // check for required tags
    if (header->version.empty())
        header->version = string(hts_version());
}

void SamHeaderCodec::DecodeSequenceLine(const string& line,
                                        SamHeader* header)
{
    static string Token_SN = string("SN");
    static string Token_LN = string("LN");
    static string Token_AS = string("AS");
    static string Token_M5 = string("M5");
    static string Token_SP = string("SP");
    static string Token_UR = string("UR");

    // split SQ line into tokens
    const vector<string>& tokens = Split(line);
    if (tokens.empty())
        return;

    SamSequence seq;

    // iterate over tokens
    auto tokenEnd = tokens.cend();
    for (auto tokenIter = tokens.cbegin(); tokenIter != tokenEnd; ++tokenIter) {

        // get tag/value
        const string tokenTag   = (*tokenIter).substr(0,2);
        const string tokenValue = (*tokenIter).substr(3);

        // set sequence contents
        if      (tokenTag == Token_SN) seq.name = tokenValue;
        else if (tokenTag == Token_LN) seq.length = tokenValue;
        else if (tokenTag == Token_AS) seq.assemblyId = tokenValue;
        else if (tokenTag == Token_M5) seq.checksum = tokenValue;
        else if (tokenTag == Token_SP) seq.species = tokenValue;
        else if (tokenTag == Token_UR) seq.uri = tokenValue;
    }

    // TODO: warn for missing values?
    if (!seq.name.empty() && !seq.length.empty())
        header->sequences.Add(seq);
}

void SamHeaderCodec::DecodeReadGroupLine(const string& line,
                                         SamHeader* header)
{
    static string Token_ID = string("ID");
    static string Token_CN = string("CN");
    static string Token_DS = string("DS");
    static string Token_DT = string("DT");
    static string Token_FO = string("FO");
    static string Token_KS = string("KS");
    static string Token_LB = string("LB");
    static string Token_PG = string("PG");
    static string Token_PI = string("PI");
    static string Token_PL = string("PL");
    static string Token_PU = string("PU");
    static string Token_SM = string("SM");

    // split SQ line into tokens
    const vector<string>& tokens = Split(line);
    if (tokens.empty())
        return;

    SamReadGroup rg;

    // iterate over tokens
    auto tokenEnd  = tokens.cend();
    for (auto tokenIter = tokens.cbegin(); tokenIter != tokenEnd; ++tokenIter) {

        // get tag/value
        const string tokenTag   = (*tokenIter).substr(0,2);
        const string tokenValue = (*tokenIter).substr(3);

        // set sequence contents
        if      (tokenTag == Token_ID) rg.id = tokenValue;
        else if (tokenTag == Token_CN) rg.sequencingCenter = tokenValue;
        else if (tokenTag == Token_DS) rg.description = tokenValue;
        else if (tokenTag == Token_DT) rg.date = tokenValue;
        else if (tokenTag == Token_FO) rg.flowOrder = tokenValue;
        else if (tokenTag == Token_KS) rg.keySequence = tokenValue;
        else if (tokenTag == Token_LB) rg.library = tokenValue;
        else if (tokenTag == Token_PG) rg.programs = tokenValue;
        else if (tokenTag == Token_PI) rg.predictedInsertSize = tokenValue;
        else if (tokenTag == Token_PL) rg.platform = tokenValue;
        else if (tokenTag == Token_PU) rg.platformUnit = tokenValue;
        else if (tokenTag == Token_SM) rg.sample = tokenValue;
    }

    // TODO: warn for missing values?
    if (!rg.id.empty())
        header->readGroups.Add(rg);
}

void SamHeaderCodec::DecodeProgramLine(const string& line,
                                       SamHeader* header)
{
    const string Token_ID = string("ID");
    const string Token_CL = string("CL");
    const string Token_DS = string("DS");
    const string Token_PN = string("PN");
    const string Token_PP = string("PP");
    const string Token_VN = string("VN");

    // split SQ line into tokens
    const vector<string>& tokens = Split(line);
    if (tokens.empty())
        return;

    SamProgram prog;

    // iterate over tokens
    auto tokenEnd  = tokens.cend();
    for (auto tokenIter = tokens.cbegin(); tokenIter != tokenEnd; ++tokenIter) {

        // get tag/value
        const string tokenTag   = (*tokenIter).substr(0,2);
        const string tokenValue = (*tokenIter).substr(3);

        // set sequence contents
        if      (tokenTag == Token_ID) prog.id = tokenValue;
        else if (tokenTag == Token_CL) prog.commandLine = tokenValue;
        else if (tokenTag == Token_DS) prog.description = tokenValue;
        else if (tokenTag == Token_PN) prog.name = tokenValue;
        else if (tokenTag == Token_PP) prog.previousProgramId = tokenValue;
        else if (tokenTag == Token_VN) prog.version = tokenValue;
    }

    // TODO: warn for missing values?
    if (!prog.id.empty())
        header->programs.Add(prog);
}

void SamHeaderCodec::DecodeCommentLine(const string& line,
                                       SamHeader* header)
{
    header->comments.push_back(line);
}

string SamHeaderCodec::Encode(const SamHeader& header)
{
    // clear out stream
    stringstream out("");

    // @HD
    EncodeHeaderLine(header.version, header.sortOrder, header.pacbioBamVersion, &out);

    // @SQ
    auto seqEnd = header.sequences.ConstEnd();
    for (auto seqIter = header.sequences.ConstBegin(); seqIter != seqEnd; ++seqIter)
        EncodeSequence(*seqIter, &out);

    // @RG
    auto rgEnd  = header.readGroups.ConstEnd();
    for (auto rgIter = header.readGroups.ConstBegin(); rgIter != rgEnd; ++rgIter)
        EncodeReadGroup(*rgIter, &out);

    // @PG
    auto progEnd = header.programs.ConstEnd();
    for (auto progIter = header.programs.ConstBegin(); progIter != progEnd; ++progIter)
        EncodeProgram(*progIter, &out);

    // @CO
    auto cEnd  = header.comments.cend();
    for (auto cIter = header.comments.cbegin(); cIter != cEnd; ++cIter)
        EncodeComment(*cIter, &out);

    // return result
    return out.str();
}

static inline
const string FormatTag(const string& tag, const string& value) {
    return string('\t' + tag + ':' + value);
}

void SamHeaderCodec::EncodeHeaderLine(const string& version,
                                      const string& sortOrder,
                                      const string& pbVersion,
                                      stringstream* out)
{
    const string& outputVersion   = (version.empty()   ? string(hts_version()) : version);
    const string& outputSortOrder = (sortOrder.empty() ? string("unknown") : sortOrder);
//    const string& outputPBVersion = (pbVersion.empty() ? string("3.0b3") : pbVersion);

    *out << "@HD"
         << FormatTag("VN", outputVersion)
         << FormatTag("SO", outputSortOrder);

    if (!pbVersion.empty())
         *out << FormatTag("pb", pbVersion);

     *out << endl;
}

void SamHeaderCodec::EncodeReadGroup(const SamReadGroup& readGroup,
                                     stringstream* out)
{
    *out << "@RG"
         << FormatTag("ID", readGroup.id);

    if (!readGroup.sequencingCenter.empty())    *out << FormatTag("CN", readGroup.sequencingCenter);
    if (!readGroup.description.empty())         *out << FormatTag("DS", readGroup.description);
    if (!readGroup.date.empty())                *out << FormatTag("DT", readGroup.date);
    if (!readGroup.flowOrder.empty())           *out << FormatTag("FL", readGroup.flowOrder);
    if (!readGroup.keySequence.empty())         *out << FormatTag("KS", readGroup.keySequence);
    if (!readGroup.library.empty())             *out << FormatTag("LB", readGroup.library);
    if (!readGroup.programs.empty())            *out << FormatTag("PG", readGroup.programs);
    if (!readGroup.predictedInsertSize.empty()) *out << FormatTag("PI", readGroup.predictedInsertSize);
    if (!readGroup.platform.empty())            *out << FormatTag("PL", readGroup.platform);
    if (!readGroup.platformUnit.empty())        *out << FormatTag("PU", readGroup.platformUnit);
    if (!readGroup.sample.empty())              *out << FormatTag("SM", readGroup.sample);

    *out << endl;
}

void SamHeaderCodec::EncodeSequence(const SamSequence& sequence,
                                    stringstream* out)
{
    *out << "@SQ"
         << FormatTag("SN", sequence.name);

    if (!sequence.length.empty())     *out << FormatTag("LN", sequence.length);
    if (!sequence.assemblyId.empty()) *out << FormatTag("AS", sequence.assemblyId);
    if (!sequence.checksum.empty())   *out << FormatTag("M5", sequence.checksum);
    if (!sequence.species.empty())    *out << FormatTag("SP", sequence.species);
    if (!sequence.uri.empty())        *out << FormatTag("UR", sequence.uri);

    *out << endl;
}

void SamHeaderCodec::EncodeProgram(const SamProgram& program,
                                   stringstream* out)
{
    *out << "@PG"
         << FormatTag("ID", program.id);

    if (!program.name.empty())              *out << FormatTag("PN", program.name);
    if (!program.version.empty())           *out << FormatTag("VN", program.version);
    if (!program.description.empty())       *out << FormatTag("DS", program.description);
    if (!program.previousProgramId.empty()) *out << FormatTag("PP", program.previousProgramId);
    if (!program.commandLine.empty())       *out << FormatTag("CL", program.commandLine);

    *out << endl;
}

void SamHeaderCodec::EncodeComment(const string& comment,
                                   stringstream* out)
{
    *out << "@CO\t" << comment << endl;
}
