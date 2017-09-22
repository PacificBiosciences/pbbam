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
/// \file ProgramInfo.h
/// \brief Defines the ProgramInfo class.
//
// Author: Derek Barnett

#ifndef PROGRAMINFO_H
#define PROGRAMINFO_H

#include "pbbam/Config.h"
#include <map>
#include <string>

namespace PacBio {
namespace BAM {

/// \brief The ProgramInfo class represents a program entry (\@PG) in the SAM
///        header.
///
class PBBAM_EXPORT ProgramInfo
{
public:
    /// \name Conversion & Validation
    ///

    /// \brief Creates a ProgramInfo object from SAM-formatted text.
    ///
    /// \param[in] sam  SAM-formatted text
    /// \returns program info object
    ///
    static ProgramInfo FromSam(const std::string& sam);

    /// \brief Converts a ProgramInfo object to its SAM-formatted text.
    ///
    /// \param[in] prog     input ProgramInfo object
    /// \returns SAM-formatted text (no trailing newline)
    ///
    static std::string ToSam(const ProgramInfo& prog);

    /// \}

public:
    /// \name Constructors & Related Methods
    /// \{

    /// \brief Creates an empty program info object.
    ProgramInfo() = default;

    /// \brief Creates a program info object with an ID.
    ///
    /// \param[in] id       program ID (\@PG:ID)
    ///
    ProgramInfo(std::string id);


    ProgramInfo(const ProgramInfo& other) = default;
    ProgramInfo(ProgramInfo&& other) = default;
    ProgramInfo& operator=(const ProgramInfo& other) = default;
    ProgramInfo& operator=(ProgramInfo&& other) = default;
    ~ProgramInfo() = default;

    /// \}

public:
    /// \name Conversion & Validation
    ///

    /// \returns true if program info is valid
    ///
    /// Currently this checks to see that ProgramInfo::Id does not contain an
    /// empty string.
    ///
    bool IsValid() const;

    /// \brief Converts this object to its SAM-formatted text.
    ///
    /// \returns SAM-formatted text (no trailing newline)
    ///
    std::string ToSam() const;

    /// \}

public:
    /// \name Attributes
    /// \{

    /// \returns string value of \@PG:CL
    std::string CommandLine() const;

    /// \returns any non-standard tags added to the \@PG entry
    ///
    /// Result map consists of {tagName => value}.
    ///
    std::map<std::string, std::string> CustomTags() const;

    /// \returns string value of \@PG:DS
    std::string Description() const;

    /// \returns string value of \@PG:ID
    std::string Id() const;

    /// \returns string value of \@PG:PN
    std::string Name() const;

    /// \returns string value of \@PG:PP
    std::string PreviousProgramId() const;

    /// \returns string value of \@PG:VN
    std::string Version() const;

    /// \}

public:
    /// \name Attributes
    /// \{

    /// \brief Sets the value for \@PG:CL
    ///
    /// \param[in] cmd      new value
    /// \returns reference to this object
    ///
    ProgramInfo& CommandLine(const std::string& cmd);

    /// \brief Sets a new collection of non-standard tags.
    ///
    /// Custom tag map entries should consist of {tagName => value}.
    ///
    /// \param[in] custom      new tags
    /// \returns reference to this object
    ///
    ProgramInfo& CustomTags(const std::map<std::string, std::string>& custom);

    /// \brief Sets the value for \@PG:DS
    ///
    /// \param[in] description      new value
    /// \returns reference to this object
    ///
    ProgramInfo& Description(const std::string& description);

    /// \brief Sets the value for \@PG:ID
    ///
    /// \param[in] id      new value
    /// \returns reference to this object
    ///
    ProgramInfo& Id(const std::string& id);

    /// \brief Sets the value for \@PG:PN
    ///
    /// \param[in] name      new value
    /// \returns reference to this object
    ///
    ProgramInfo& Name(const std::string& name);

    /// \brief Sets the value for \@PG:PP
    ///
    /// \param[in] id      new value
    /// \returns reference to this object
    ///
    ProgramInfo& PreviousProgramId(const std::string& id);

    /// \brief Sets the value for \@PG:VN
    ///
    /// \param[in] version      new value
    /// \returns reference to this object
    ///
    ProgramInfo& Version(const std::string& version);

    /// \}

private:
    std::string commandLine_;       // CL:<CommandLine>
    std::string description_;       // DS:<Description>
    std::string id_;                // ID:<ID>  * must be unique for valid SAM *
    std::string name_;              // PN:<Name>
    std::string previousProgramId_; // PP:<PreviousProgramID>
    std::string version_;           // VN:<Version>

    // custom attributes
    std::map<std::string, std::string> custom_;     // tag => value
};

} // namespace BAM
} // namespace PacBio

#include "pbbam/internal/ProgramInfo.inl"

#endif // PROGRAMINFO_H
