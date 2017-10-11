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
/// \file SequenceInfo.inl
/// \brief Inline implementations for the SequenceInfo class.
//
// Author: Derek Barnett

#include "pbbam/SequenceInfo.h"

namespace PacBio {
namespace BAM {

inline bool SequenceInfo::operator==(const SequenceInfo& other) const
{
    return assemblyId_ == other.assemblyId_ &&
           checksum_   == other.checksum_   &&
           length_     == other.length_     &&
           name_       == other.name_       &&
           species_    == other.species_    &&
           uri_        == other.uri_        &&
           custom_     == other.custom_;
}

inline bool SequenceInfo::operator!=(const SequenceInfo& other) const
{ return !(*this == other); }

inline std::string SequenceInfo::AssemblyId() const
{ return assemblyId_; }

inline SequenceInfo& SequenceInfo::AssemblyId(const std::string& id)
{ assemblyId_ = id; return *this; }

inline std::string SequenceInfo::Checksum() const
{ return checksum_; }

inline SequenceInfo& SequenceInfo::Checksum(const std::string& checksum)
{ checksum_ = checksum; return *this; }

inline std::map<std::string, std::string> SequenceInfo::CustomTags() const
{ return custom_; }

inline SequenceInfo& SequenceInfo::CustomTags(const std::map<std::string, std::string>& custom)
{ custom_ = custom; return *this; }

inline std::string SequenceInfo::Length() const
{ return length_; }

inline SequenceInfo& SequenceInfo::Length(const std::string& length)
{ length_ = length; return *this; }

inline std::string SequenceInfo::Name() const
{ return name_; }

inline SequenceInfo& SequenceInfo::Name(const std::string& name)
{ name_ = name; return *this; }

inline std::string SequenceInfo::Species() const
{ return species_; }

inline SequenceInfo& SequenceInfo::Species(const std::string& species)
{ species_ = species; return *this; }

inline std::string SequenceInfo::ToSam(const SequenceInfo& seq)
{ return seq.ToSam(); }

inline std::string SequenceInfo::Uri() const
{ return uri_; }

inline SequenceInfo& SequenceInfo::Uri(const std::string& uri)
{ uri_ = uri; return *this; }

} // namespace BAM
} // namespace PacBio
