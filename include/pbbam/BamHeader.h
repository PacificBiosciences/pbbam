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

#ifndef BAMHEADER_H
#define BAMHEADER_H

#include "pbbam/Config.h"
#include "pbbam/ProgramInfo.h"
#include "pbbam/ReadGroupInfo.h"
#include "pbbam/SequenceInfo.h"
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

class PBBAM_EXPORT BamHeader
{
public:
    typedef PBBAM_SHARED_PTR<BamHeader> SharedPtr;

public:
    /// \name Conversion Methods
    /// \{

    static BamHeader::SharedPtr FromSam(const std::string& sam);

    /// \}

public:
    /// \name Constructors & Related Methods
    /// \{

    BamHeader(void);
    BamHeader(const std::string& text);
    BamHeader(const BamHeader& other);
    BamHeader(BamHeader&& other);
    BamHeader& operator=(const BamHeader& other);
    BamHeader& operator=(BamHeader&& other);
    ~BamHeader(void);

    /// \}

public:
    /// \name General
    /// \{

    std::string PacBioBamVersion(void) const;
    std::string SortOrder(void) const;
    std::string Version(void) const;

    /// \}

    /// \name Read Groups
    /// \{

    bool HasReadGroup(const std::string& id) const;
    ReadGroupInfo ReadGroup(const std::string& id) const;
    std::vector<std::string> ReadGroupIds(void) const;
    std::vector<ReadGroupInfo> ReadGroups(void) const;

    /// \}

    /// \name Sequences
    /// \{

    bool HasSequence(const std::string& name) const;
    int32_t SequenceId(const std::string& name) const;
    std::string SequenceLength(const int32_t id) const;
    std::string SequenceName(const int32_t id) const;
    std::vector<std::string> SequenceNames(void) const;
    SequenceInfo Sequence(const std::string& name) const;
    std::vector<SequenceInfo> Sequences(void) const;

    /// \}

    /// \name Programs
    /// \{

    bool HasProgram(const std::string& id) const;
    ProgramInfo Program(const std::string& id) const;
    std::vector<std::string> ProgramIds(void) const;
    std::vector<ProgramInfo> Programs(void) const;

    /// \}

    /// \name Comments
    /// \{

    std::vector<std::string> Comments(void) const;

    /// \}

    /// \name Conversion Methods
    /// \{

    std::string ToSam(void) const;

    /// \}

public:

    /// \name General
    /// \{

    BamHeader& PacBioBamVersion(const std::string& version);
    BamHeader& SortOrder(const std::string& order);
    BamHeader& Version(const std::string& version);

    /// \}

    /// \name Read Groups
    /// \{

    BamHeader& AddReadGroup(const ReadGroupInfo& readGroup);
    BamHeader& ClearReadGroups(void);
    BamHeader& ReadGroups(const std::vector<ReadGroupInfo>& readGroups);

    /// \}

    /// \name Sequences
    /// \{

    BamHeader& AddSequence(const SequenceInfo& sequence);
    BamHeader& ClearSequences(void);
    BamHeader& Sequences(const std::vector<SequenceInfo>& sequences);

    /// \}

    /// \name Programs
    /// \{

    BamHeader& AddProgram(const ProgramInfo& pg);
    BamHeader& ClearPrograms(void);
    BamHeader& Programs(const std::vector<ProgramInfo>& programs);

    /// \}

    /// \name Comments
    /// \{

    BamHeader& AddComment(const std::string& comment);
    BamHeader& ClearComments(void);
    BamHeader& Comments(const std::vector<std::string>& comments);

    /// \}

private:
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

inline BamHeader& BamHeader::AddComment(const std::string& comment)
{ comments_.push_back(comment); return *this; }

inline BamHeader& BamHeader::AddProgram(const ProgramInfo& pg)
{ programs_[pg.Id()] = pg; return *this; }

inline BamHeader& BamHeader::AddReadGroup(const ReadGroupInfo& readGroup)
{ readGroups_[readGroup.Id()] = readGroup; return *this; }

inline BamHeader& BamHeader::ClearComments(void)
{ comments_.clear(); return* this; }

inline BamHeader& BamHeader::ClearPrograms(void)
{ programs_.clear(); return *this; }

inline BamHeader& BamHeader::ClearReadGroups(void)
{ readGroups_.clear(); return *this; }

inline BamHeader& BamHeader::ClearSequences(void)
{
    sequenceIdLookup_.clear();
    sequences_.clear();
    return *this;
}

inline std::vector<std::string> BamHeader::Comments(void) const
{ return comments_; }

inline BamHeader& BamHeader::Comments(const std::vector<std::string>& comments)
{ comments_ = comments; return *this; }

inline BamHeader::SharedPtr BamHeader::FromSam(const std::string& sam)
{ return BamHeader::SharedPtr(new BamHeader(sam)); }

inline bool BamHeader::HasProgram(const std::string& id) const
{ return programs_.find(id) != programs_.cend(); }

inline bool BamHeader::HasReadGroup(const std::string& id) const
{ return readGroups_.find(id) != readGroups_.cend(); }

inline bool BamHeader::HasSequence(const std::string& name) const
{ return sequenceIdLookup_.find(name) != sequenceIdLookup_.cend(); }

inline std::string BamHeader::PacBioBamVersion(void) const
{ return pacbioBamVersion_; }

inline BamHeader& BamHeader::PacBioBamVersion(const std::string& version)
{ pacbioBamVersion_ = version; return *this; }

inline ProgramInfo BamHeader::Program(const std::string& id) const
{
    const auto iter = programs_.find(id);
    if (iter == programs_.cend())
        return ProgramInfo();
    return iter->second;
}

inline ReadGroupInfo BamHeader::ReadGroup(const std::string& id) const
{
    const auto iter = readGroups_.find(id);
    if (iter == readGroups_.cend())
        return ReadGroupInfo();
    return iter->second;
}

inline int32_t BamHeader::SequenceId(const std::string& name) const
{
    const auto iter = sequenceIdLookup_.find(name);
    if (iter == sequenceIdLookup_.cend())
        return -1;
    return iter->second;
}

inline std::string BamHeader::SequenceLength(const int32_t id) const
{
    if (id < 0 || (size_t)id >= sequences_.size())
        return 0;
    return sequences_.at(id).Length();
}

inline std::string BamHeader::SequenceName(const int32_t id) const
{
    if (id < 0 || (size_t)id >= sequences_.size())
        return std::string();
    return sequences_.at(id).Name();
}

inline std::vector<SequenceInfo> BamHeader::Sequences(void) const
{ return sequences_; }

inline std::string BamHeader::SortOrder(void) const
{ return sortOrder_; }

inline BamHeader& BamHeader::SortOrder(const std::string& order)
{ sortOrder_ = order; return *this; }

inline std::string BamHeader::Version(void) const
{ return version_; }

inline BamHeader& BamHeader::Version(const std::string& version)
{ version_ = version; return *this; }

} // namespace BAM
} // namespace PacBio

#endif // BAMHEADER_H
