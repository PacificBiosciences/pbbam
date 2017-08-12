
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

namespace PacBio {
namespace BAM {
namespace internal {

static std::string ProgramInfoTokenID{"ID"};
static std::string ProgramInfoTokenCL{"CL"};
static std::string ProgramInfoTokenDS{"DS"};
static std::string ProgramInfoTokenPN{"PN"};
static std::string ProgramInfoTokenPP{"PP"};
static std::string ProgramInfoTokenVN{"VN"};

} // namespace internal

ProgramInfo::ProgramInfo(const std::string& id)
    : id_(id)
{ }

ProgramInfo ProgramInfo::FromSam(const std::string& sam)
{
    // pop off '@PG\t', then split rest of line into tokens
    const std::vector<std::string>& tokens = internal::Split(sam.substr(4), '\t');
    if (tokens.empty())
        return ProgramInfo();

    ProgramInfo prog;
    std::map<std::string, std::string> custom;

    // iterate over tokens
    for (const std::string& token : tokens) {
        const std::string& tokenTag   = token.substr(0,2);
        const std::string& tokenValue = token.substr(3);

        // set program contents
        if      (tokenTag == internal::ProgramInfoTokenID) prog.Id(tokenValue);
        else if (tokenTag == internal::ProgramInfoTokenCL) prog.CommandLine(tokenValue);
        else if (tokenTag == internal::ProgramInfoTokenDS) prog.Description(tokenValue);
        else if (tokenTag == internal::ProgramInfoTokenPN) prog.Name(tokenValue);
        else if (tokenTag == internal::ProgramInfoTokenPP) prog.PreviousProgramId(tokenValue);
        else if (tokenTag == internal::ProgramInfoTokenVN) prog.Version(tokenValue);

        // otherwise, "custom" tag
        else
            custom[tokenTag] = tokenValue;
    }

    prog.CustomTags(custom);
    return prog;
}

std::string ProgramInfo::ToSam(void) const
{
    std::stringstream out;
    out << "@PG"
        << internal::MakeSamTag(internal::ProgramInfoTokenID, id_);

    if (!name_.empty())              out << internal::MakeSamTag(internal::ProgramInfoTokenPN, name_);
    if (!version_.empty())           out << internal::MakeSamTag(internal::ProgramInfoTokenVN, version_);
    if (!description_.empty())       out << internal::MakeSamTag(internal::ProgramInfoTokenDS, description_);
    if (!previousProgramId_.empty()) out << internal::MakeSamTag(internal::ProgramInfoTokenPP, previousProgramId_);
    if (!commandLine_.empty())       out << internal::MakeSamTag(internal::ProgramInfoTokenCL, commandLine_);

    // append any custom tags
    std::map<std::string, std::string>::const_iterator customIter = custom_.cbegin();
    std::map<std::string, std::string>::const_iterator customEnd  = custom_.cend();
    for ( ; customIter != customEnd; ++customIter )
        out << internal::MakeSamTag(customIter->first, customIter->second);

    return out.str();
}

} // namespace BAM
} // namespace PacBio
