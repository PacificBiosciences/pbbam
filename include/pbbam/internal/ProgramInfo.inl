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
/// \file ProgramInfo.inl
/// \brief Inline implementations for the ProgramInfo class.
//
// Author: Derek Barnett

#include "pbbam/ProgramInfo.h"

namespace PacBio {
namespace BAM {

inline std::string ProgramInfo::CommandLine() const
{ return commandLine_; }

inline ProgramInfo& ProgramInfo::CommandLine(const std::string& cmd)
{ commandLine_ = cmd; return *this; }

inline std::map<std::string, std::string> ProgramInfo::CustomTags() const
{ return custom_; }

inline ProgramInfo& ProgramInfo::CustomTags(const std::map<std::string,
                                            std::string>& custom)
{ custom_ = custom; return *this; }

inline std::string ProgramInfo::Description() const
{ return description_; }

inline ProgramInfo& ProgramInfo::Description(const std::string& description)
{ description_ = description; return *this; }

inline std::string ProgramInfo::Id() const
{ return id_; }

inline ProgramInfo& ProgramInfo::Id(const std::string& id)
{ id_ = id; return *this; }

inline bool ProgramInfo::IsValid() const
{ return !id_.empty(); }

inline std::string ProgramInfo::Name() const
{ return name_; }

inline ProgramInfo& ProgramInfo::Name(const std::string& name)
{ name_ = name; return *this; }

inline std::string ProgramInfo::PreviousProgramId() const
{ return previousProgramId_; }

inline ProgramInfo& ProgramInfo::PreviousProgramId(const std::string& id)
{ previousProgramId_ = id; return *this; }

inline std::string ProgramInfo::ToSam(const ProgramInfo& prog)
{ return prog.ToSam(); }

inline std::string ProgramInfo::Version() const
{ return version_; }

inline ProgramInfo& ProgramInfo::Version(const std::string& version)
{ version_ = version; return *this; }

} // namespace BAM
} // namespace PacBio
