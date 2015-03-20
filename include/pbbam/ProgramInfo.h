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

// Author: Derek Barnett

#ifndef PROGRAMINFO_H
#define PROGRAMINFO_H

#include "pbbam/Config.h"
#include <map>
#include <string>

namespace PacBio {
namespace BAM {

class PBBAM_EXPORT ProgramInfo
{
public:
    /// \name Conversion & Validation
    ///

    static ProgramInfo FromSam(const std::string& sam);

    static std::string ToSam(const ProgramInfo& prog);

    /// \}

public:
    /// \name Constructors & Related Methods
    /// \{

    ProgramInfo(void);
    ProgramInfo(const std::string& id);
    ProgramInfo(const ProgramInfo& other);
    ProgramInfo(ProgramInfo&& other);
    ProgramInfo& operator=(const ProgramInfo& other);
    ProgramInfo& operator=(ProgramInfo&& other);
    ~ProgramInfo(void);

    /// \}

public:
    /// \name Attributes
    /// \{

    std::string CommandLine(void) const;

    std::map<std::string, std::string> CustomTags(void) const;

    std::string Description(void) const;

    std::string Id(void) const;

    std::string Name(void) const;

    std::string PreviousProgramId(void) const;

    std::string Version(void) const;

    /// \}

    /// \name Conversion & Validation
    ///

    bool IsValid(void) const;

    std::string ToSam(void) const;

    /// \}

public:
    /// \name Attributes
    /// \{

    ProgramInfo& CommandLine(const std::string& cmd);

    ProgramInfo& CustomTags(const std::map<std::string, std::string>& custom);

    ProgramInfo& Description(const std::string& description);

    ProgramInfo& Id(const std::string& id);

    ProgramInfo& Name(const std::string& name);

    ProgramInfo& PreviousProgramId(const std::string& id);

    ProgramInfo& Version(const std::string& version);

    /// \}

private:
    std::string commandLine_;            // CL:<CommandLine>
    std::string description_;            // DS:<Description>
    std::string id_;                     // ID:<ID>              * Unique ID required for valid SAM header*
    std::string name_;                   // PN:<Name>
    std::string previousProgramId_;      // PP:<PreviousProgramID>
    std::string version_;                // VN:<Version>

    // custom attributes
    std::map<std::string, std::string> custom_; // tag => value
};

inline std::string ProgramInfo::CommandLine(void) const
{ return commandLine_; }

inline ProgramInfo& ProgramInfo::CommandLine(const std::string& cmd)
{ commandLine_ = cmd; return *this; }

inline std::map<std::string, std::string> ProgramInfo::CustomTags(void) const
{ return custom_; }

inline ProgramInfo& ProgramInfo::CustomTags(const std::map<std::string, std::string>& custom)
{ custom_ = custom; return *this; }

inline std::string ProgramInfo::Description(void) const
{ return description_; }

inline ProgramInfo& ProgramInfo::Description(const std::string& description)
{ description_ = description; return *this; }

inline std::string ProgramInfo::Id(void) const
{ return id_; }

inline ProgramInfo& ProgramInfo::Id(const std::string& id)
{ id_ = id; return *this; }

inline bool ProgramInfo::IsValid(void) const
{ return !id_.empty(); }

inline std::string ProgramInfo::Name(void) const
{ return name_; }

inline ProgramInfo& ProgramInfo::Name(const std::string& name)
{ name_ = name; return *this; }

inline std::string ProgramInfo::PreviousProgramId(void) const
{ return previousProgramId_; }

inline ProgramInfo& ProgramInfo::PreviousProgramId(const std::string& id)
{ previousProgramId_ = id; return *this; }

inline std::string ProgramInfo::ToSam(const ProgramInfo& prog)
{ return prog.ToSam(); }

inline std::string ProgramInfo::Version(void) const
{ return version_; }

inline ProgramInfo& ProgramInfo::Version(const std::string& version)
{ version_ = version; return *this; }

} // namespace BAM
} // namespace PacBio

#endif // PROGRAMINFO_H
