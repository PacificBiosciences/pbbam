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
/// \file BamHeader.inl
/// \brief Inline implementations for the BamHeader class.
//
// Author: Derek Barnett

#include "pbbam/BamHeader.h"

namespace PacBio {
namespace BAM {
namespace internal {

class BamHeaderPrivate
{
public:
    std::string version_;
    std::string pacbioBamVersion_;
    std::string sortOrder_;
    std::map<std::string, std::string> headerLineCustom_;

    std::map<std::string, ReadGroupInfo> readGroups_; // id => read group info
    std::map<std::string, ProgramInfo> programs_;     // id => program info
    std::vector<std::string> comments_;

    // we need to preserve insertion order, use lookup for access by name
    std::vector<SequenceInfo> sequences_;
    std::map<std::string, int32_t> sequenceIdLookup_;
};

} // namespace internal

inline BamHeader::BamHeader()
    : d_(new internal::BamHeaderPrivate)
{ }

inline BamHeader BamHeader::operator+(const BamHeader& other) const
{ return DeepCopy() += other; }

inline BamHeader& BamHeader::AddComment(const std::string& comment)
{ d_->comments_.push_back(comment); return *this; }

inline BamHeader& BamHeader::AddProgram(const ProgramInfo& pg)
{ d_->programs_[pg.Id()] = pg; return *this; }

inline BamHeader& BamHeader::AddReadGroup(const ReadGroupInfo& readGroup)
{ d_->readGroups_[readGroup.Id()] = readGroup; return *this; }

inline BamHeader& BamHeader::ClearComments()
{ d_->comments_.clear(); return* this; }

inline BamHeader& BamHeader::ClearPrograms()
{ d_->programs_.clear(); return *this; }

inline BamHeader& BamHeader::ClearReadGroups()
{ d_->readGroups_.clear(); return *this; }

inline std::vector<std::string> BamHeader::Comments() const
{ return d_->comments_; }

inline BamHeader& BamHeader::Comments(const std::vector<std::string>& comments)
{ d_->comments_ = comments; return *this; }

inline bool BamHeader::HasProgram(const std::string& id) const
{ return d_->programs_.find(id) != d_->programs_.cend(); }

inline bool BamHeader::HasReadGroup(const std::string& id) const
{ return d_->readGroups_.find(id) != d_->readGroups_.cend(); }

inline bool BamHeader::HasSequence(const std::string& name) const
{ return d_->sequenceIdLookup_.find(name) != d_->sequenceIdLookup_.cend(); }

inline size_t BamHeader::NumSequences() const
{ return d_->sequences_.size(); }

inline std::string BamHeader::PacBioBamVersion() const
{ return d_->pacbioBamVersion_; }

inline SequenceInfo BamHeader::Sequence(const int32_t id) const
{ return d_->sequences_.at(id); }

inline std::string BamHeader::SequenceLength(const int32_t id) const
{ return Sequence(id).Length(); }

inline std::string BamHeader::SequenceName(const int32_t id) const
{ return Sequence(id).Name(); }

inline std::vector<SequenceInfo> BamHeader::Sequences() const
{ return d_->sequences_; }

inline std::string BamHeader::SortOrder() const
{ return d_->sortOrder_; }

inline BamHeader& BamHeader::SortOrder(const std::string& order)
{ d_->sortOrder_ = order; return *this; }

inline std::string BamHeader::Version() const
{ return d_->version_; }

inline BamHeader& BamHeader::Version(const std::string& version)
{ d_->version_ = version; return *this; }

} // namespace BAM
} // namespace PacBio
