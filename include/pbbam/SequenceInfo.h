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
/// \file SequenceInfo.h
/// \brief Defines the SequenceInfo class.
//
// Author: Derek Barnett

#ifndef SEQUENCEINFO_H
#define SEQUENCEINFO_H

#include "pbbam/Config.h"
#include <map>
#include <string>

namespace PacBio {
namespace BAM {

/// \brief The SequenceInfo class represents a program entry (\@SQ) in the SAM
///        header.
///
class PBBAM_EXPORT SequenceInfo
{
public:
    /// \name Conversion & Validation
    ///

    /// \brief Creates a SequenceInfo object from SAM-formatted text.
    ///
    /// \param[in] sam  SAM-formatted text
    /// \returns program info object
    ///
    static SequenceInfo FromSam(const std::string& sam);

    /// \brief Converts a SequenceInfo object to its SAM-formatted text.
    ///
    /// \param[in] seq     input SequenceInfo object
    /// \returns SAM-formatted text (no trailing newline)
    ///
    static std::string ToSam(const SequenceInfo& seq);

    /// \}

public:
    /// \name Constructors & Related Methods
    /// \{

    /// \brief Creates a sequence info object with name & (optional) length.
    ///
    /// \param[in] name       sequence name (\@SQ:SN)
    /// \param[in] length     sequence length (\@SQ:LN)
    ///
    SequenceInfo(std::string name, std::string length = "0");

    SequenceInfo() = default;
    SequenceInfo(const SequenceInfo&) = default;
    SequenceInfo(SequenceInfo&&) = default;
    SequenceInfo& operator=(const SequenceInfo&) = default;
    SequenceInfo& operator=(SequenceInfo&&) = default;
    ~SequenceInfo() = default;

    /// \}

public:
    /// \name Operators
    /// \{

    bool operator==(const SequenceInfo& other) const;
    bool operator!=(const SequenceInfo& other) const;

    /// \}

public:
    /// \name Conversion & Validation
    ///

    /// \returns true if sequence info is valid
    ///
    /// Currently this checks to see that Name is non-empty and Length is within
    /// the accepted range.
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

    /// \returns string value of \@SQ:AS
    std::string AssemblyId() const;

    /// \returns string value of \@SQ:M5
    std::string Checksum() const;

    /// \returns any non-standard tags added to the \@PG entry
    ///
    /// Result map consists of {tagName => value}.
    ///
    std::map<std::string, std::string> CustomTags() const;

    /// \returns string value of \@SQ:LN
    std::string Length() const;

    /// \returns string value of \@SQ:SN
    std::string Name() const;

    /// \returns string value of \@SQ:SP
    std::string Species() const;

    /// \returns string value of \@SQ:UR
    std::string Uri() const;

    /// \}

public:
    /// \name Attributes
    /// \{

    /// \brief Sets the value for \@SQ:AS
    ///
    /// \param[in] id      new value
    /// \returns reference to this object
    ///
    SequenceInfo& AssemblyId(const std::string& id);

    /// \brief Sets the value for \@SQ:M5
    ///
    /// \param[in] checksum      new value
    /// \returns reference to this object
    ///
    SequenceInfo& Checksum(const std::string& checksum);

    /// \brief Sets a new collection of non-standard tags.
    ///
    /// Custom tag map entries should consist of {tagName => value}.
    ///
    /// \param[in] custom      new tags
    /// \returns reference to this object
    ///
    SequenceInfo& CustomTags(const std::map<std::string, std::string>& custom);

    /// \brief Sets the value for \@SQ:LN
    ///
    /// \param[in] length      new value
    /// \returns reference to this object
    ///
    SequenceInfo& Length(const std::string& length);

    /// \brief Sets the value for \@SQ:SN
    ///
    /// \param[in] name      new value
    /// \returns reference to this object
    ///
    SequenceInfo& Name(const std::string& name);

    /// \brief Sets the value for \@SQ:SP
    ///
    /// \param[in] species     new value
    /// \returns reference to this object
    ///
    SequenceInfo& Species(const std::string& species);

    /// \brief Sets the value for \@SQ:UR
    ///
    /// \param[in] uri      new value
    /// \returns reference to this object
    ///
    SequenceInfo& Uri(const std::string& uri);

    /// \}

private:
    std::string name_;          // SN:<Name>    * must be unique for valid SAM *
    std::string length_;        // LN:<Length>  * must be within [0 - 2^31-1] *
    std::string assemblyId_;    // AS:<AssemblyId>
    std::string checksum_;      // M5:<Checksum>
    std::string species_;       // SP:<Species>
    std::string uri_;           // UR:<URI>

    // custom attributes
    std::map<std::string, std::string> custom_; // tag => value
};

} // namespace BAM
} // namespace PacBio

#include "pbbam/internal/SequenceInfo.inl"

#endif // SEQUENCEINFO_H
