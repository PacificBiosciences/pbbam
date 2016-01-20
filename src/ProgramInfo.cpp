
// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
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
//
// File Description
/// \file ProgramInfo.cpp
/// \brief Implements the ProgramInfo class.
//
// Author: Derek Barnett

#include "pbbam/ProgramInfo.h"
#include "SequenceUtils.h"
#include <sstream>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

static string token_ID = string("ID");
static string token_CL = string("CL");
static string token_DS = string("DS");
static string token_PN = string("PN");
static string token_PP = string("PP");
static string token_VN = string("VN");

} // namespace internal
} // namespace BAM
} // namespace PacBio

ProgramInfo::ProgramInfo(void) { }

ProgramInfo::ProgramInfo(const std::string& id)
    : id_(id)
{ }

ProgramInfo::ProgramInfo(const ProgramInfo& other)
    : commandLine_(other.commandLine_)
    , description_(other.description_)
    , id_(other.id_)
    , name_(other.name_)
    , previousProgramId_(other.previousProgramId_)
    , version_(other.version_)
{  }

ProgramInfo::ProgramInfo(ProgramInfo&& other)
    : commandLine_(std::move(other.commandLine_))
    , description_(std::move(other.description_))
    , id_(std::move(other.id_))
    , name_(std::move(other.name_))
    , previousProgramId_(std::move(other.previousProgramId_))
    , version_(std::move(other.version_))
{ }

ProgramInfo::~ProgramInfo(void) { }

ProgramInfo& ProgramInfo::operator=(const ProgramInfo& other)
{
    commandLine_ = other.commandLine_;
    description_ = other.description_;
    id_ = other.id_;
    name_ = other.name_;
    previousProgramId_ = other.previousProgramId_;
    version_ = other.version_;
    return *this;
}

ProgramInfo& ProgramInfo::operator=(ProgramInfo&& other)
{
    commandLine_ = std::move(other.commandLine_);
    description_ = std::move(other.description_);
    id_ = std::move(other.id_);
    name_ = std::move(other.name_);
    previousProgramId_ = std::move(other.previousProgramId_);
    version_ = std::move(other.version_);
    return *this;
}

ProgramInfo ProgramInfo::FromSam(const string& sam)
{
    // pop off '@PG\t', then split rest of line into tokens
    const vector<string>& tokens = internal::Split(sam.substr(4), '\t');
    if (tokens.empty())
        return ProgramInfo();

    ProgramInfo prog;
    map<string, string> custom;

    // iterate over tokens
    for (const string& token : tokens) {
        const string& tokenTag   = token.substr(0,2);
        const string& tokenValue = token.substr(3);

        // set program contents
        if      (tokenTag == internal::token_ID) prog.Id(tokenValue);
        else if (tokenTag == internal::token_CL) prog.CommandLine(tokenValue);
        else if (tokenTag == internal::token_DS) prog.Description(tokenValue);
        else if (tokenTag == internal::token_PN) prog.Name(tokenValue);
        else if (tokenTag == internal::token_PP) prog.PreviousProgramId(tokenValue);
        else if (tokenTag == internal::token_VN) prog.Version(tokenValue);

        // otherwise, "custom" tag
        else
            custom[tokenTag] = tokenValue;
    }

    prog.CustomTags(custom);
    return prog;
}

string ProgramInfo::ToSam(void) const
{
    stringstream out;
    out << "@PG"
        << internal::MakeSamTag(internal::token_ID, id_);

    if (!name_.empty())              out << internal::MakeSamTag(internal::token_PN, name_);
    if (!version_.empty())           out << internal::MakeSamTag(internal::token_VN, version_);
    if (!description_.empty())       out << internal::MakeSamTag(internal::token_DS, description_);
    if (!previousProgramId_.empty()) out << internal::MakeSamTag(internal::token_PP, previousProgramId_);
    if (!commandLine_.empty())       out << internal::MakeSamTag(internal::token_CL, commandLine_);

    // append any custom tags
    map<string, string>::const_iterator customIter = custom_.cbegin();
    map<string, string>::const_iterator customEnd  = custom_.cend();
    for ( ; customIter != customEnd; ++customIter )
        out << internal::MakeSamTag(customIter->first, customIter->second);

    return out.str();
}

