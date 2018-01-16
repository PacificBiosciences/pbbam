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
/// \file BamHeader.cpp
/// \brief Implements the BamHeader class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/BamHeader.h"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <set>
#include <sstream>

#include <htslib/hts.h>

#include "StringUtils.h"
#include "Version.h"

namespace PacBio {
namespace BAM {
namespace internal {

static const std::string BamHeaderPrefixHD{"@HD"};
static const std::string BamHeaderPrefixSQ{"@SQ"};
static const std::string BamHeaderPrefixRG{"@RG"};
static const std::string BamHeaderPrefixPG{"@PG"};
static const std::string BamHeaderPrefixCO{"@CO"};

static const std::string BamHeaderTokenVN{"VN"};
static const std::string BamHeaderTokenSO{"SO"};
static const std::string BamHeaderTokenpb{"pb"};

static inline bool CheckSortOrder(const std::string& lhs, const std::string& rhs)
{
    return lhs == rhs;
}

static inline bool CheckPbVersion(const std::string& lhs, const std::string& rhs)
{
    return (Version{lhs} >= Version::Minimum && Version{rhs} >= Version::Minimum);
}

static inline bool CheckSequences(const std::string& sortOrder,
                                  const std::vector<SequenceInfo>& lhs,
                                  const std::vector<SequenceInfo>& rhs)
{
    return ((sortOrder == "coordinate") ? lhs == rhs : true);
}

static void EnsureCanMerge(const BamHeader& lhs, const BamHeader& rhs)
{
    // check compatibility
    const auto sortOrderOk = CheckSortOrder(lhs.SortOrder(), rhs.SortOrder());
    const auto pbVersionOk = CheckPbVersion(lhs.PacBioBamVersion(), rhs.PacBioBamVersion());
    const auto sequencesOk = CheckSequences(lhs.SortOrder(), lhs.Sequences(), rhs.Sequences());
    if (sortOrderOk && pbVersionOk && sequencesOk) return;

    // if any checks failed, format error message & throw
    std::stringstream e;
    e << "could not merge BAM headers:" << std::endl;

    if (!sortOrderOk) {
        e << "  mismatched sort orders (@HD:SO) : (" << lhs.SortOrder() << ", " << rhs.SortOrder()
          << ")" << std::endl;
    }

    if (!pbVersionOk) {
        e << "  incompatible PacBio BAM versions (@HD:pb) : (" << lhs.PacBioBamVersion() << ", "
          << rhs.PacBioBamVersion() << ")" << std::endl;
    }

    if (!sequencesOk) e << "  mismatched sequence lists (@SQ entries)" << std::endl;

    throw std::runtime_error(e.str());
}

}  // namespace internal

BamHeader::BamHeader(const std::string& samHeaderText) : d_(new internal::BamHeaderPrivate)
{
    std::stringstream s(samHeaderText);
    std::string line("");
    std::string firstToken;
    while (getline(s, line)) {

        // skip if line is not long enough to contain true values
        if (line.length() < 5) continue;

        // determine token at beginning of line
        firstToken = line.substr(0, 3);

        if (firstToken == internal::BamHeaderPrefixHD) {

            // pop off '@HD\t', then split HD lines into tokens
            const auto tokens = internal::Split(line.substr(4), '\t');
            for (const auto& token : tokens) {
                const auto tokenTag = token.substr(0, 2);
                const auto tokenValue = token.substr(3);

                // set header contents
                if (tokenTag == internal::BamHeaderTokenVN)
                    Version(tokenValue);
                else if (tokenTag == internal::BamHeaderTokenSO)
                    SortOrder(tokenValue);
                else if (tokenTag == internal::BamHeaderTokenpb)
                    PacBioBamVersion(tokenValue);
            }

            // check for required tags
            if (Version().empty()) Version(std::string(hts_version()));
        }

        else if (firstToken == internal::BamHeaderPrefixSQ)
            AddSequence(SequenceInfo::FromSam(line));

        else if (firstToken == internal::BamHeaderPrefixRG)
            AddReadGroup(ReadGroupInfo::FromSam(line));

        else if (firstToken == internal::BamHeaderPrefixPG)
            AddProgram(ProgramInfo::FromSam(line));

        else if (firstToken == internal::BamHeaderPrefixCO)
            AddComment(line.substr(4));
    }
}

BamHeader& BamHeader::operator+=(const BamHeader& other)
{
    internal::EnsureCanMerge(*this, other);

    // merge read groups
    for (const auto& rg : other.ReadGroups()) {
        if (!HasReadGroup(rg.Id())) AddReadGroup(rg);
    }

    // merge programs
    for (const auto& pg : other.Programs()) {
        if (!HasProgram(pg.Id())) AddProgram(pg);
    }

    // merge comments
    for (const auto& comment : other.Comments())
        AddComment(comment);

    return *this;
}

BamHeader& BamHeader::AddSequence(const SequenceInfo& sequence)
{
    d_->sequences_.push_back(sequence);
    d_->sequenceIdLookup_[sequence.Name()] = d_->sequences_.size() - 1;
    return *this;
}

BamHeader& BamHeader::ClearSequences()
{
    d_->sequenceIdLookup_.clear();
    d_->sequences_.clear();
    return *this;
}

BamHeader BamHeader::DeepCopy() const
{
    BamHeader result;
    result.d_->version_ = d_->version_;
    result.d_->pacbioBamVersion_ = d_->pacbioBamVersion_;
    result.d_->sortOrder_ = d_->sortOrder_;
    result.d_->headerLineCustom_ = d_->headerLineCustom_;
    result.d_->readGroups_ = d_->readGroups_;
    result.d_->programs_ = d_->programs_;
    result.d_->comments_ = d_->comments_;
    result.d_->sequences_ = d_->sequences_;
    result.d_->sequenceIdLookup_ = d_->sequenceIdLookup_;
    return result;
}

BamHeader& BamHeader::PacBioBamVersion(const std::string& version)
{
    d_->pacbioBamVersion_ = version;
    const auto fileVersion = internal::Version{version};
    if (fileVersion < internal::Version::Minimum) {
        auto msg = std::string{"invalid PacBio BAM version number"};
        msg += ("(" + fileVersion.ToString() + ")");
        msg += std::string{"is older than the minimum supported version"};
        msg += ("(" + internal::Version::Minimum.ToString() + ")");
        throw std::runtime_error(msg);
    }
    return *this;
}

ProgramInfo BamHeader::Program(const std::string& id) const
{
    const auto iter = d_->programs_.find(id);
    if (iter == d_->programs_.cend()) throw std::runtime_error("program ID not found");
    return iter->second;
}

std::vector<std::string> BamHeader::ProgramIds() const
{
    std::vector<std::string> result;
    result.reserve(d_->programs_.size());
    const auto end = d_->programs_.cend();
    auto iter = d_->programs_.cbegin();
    for (; iter != end; ++iter)
        result.push_back(iter->first);
    return result;
}

std::vector<ProgramInfo> BamHeader::Programs() const
{
    std::vector<ProgramInfo> result;
    result.reserve(d_->programs_.size());
    const auto end = d_->programs_.cend();
    auto iter = d_->programs_.cbegin();
    for (; iter != end; ++iter)
        result.push_back(iter->second);
    return result;
}

BamHeader& BamHeader::Programs(const std::vector<ProgramInfo>& programs)
{
    d_->programs_.clear();
    for (const auto& pg : programs)
        d_->programs_[pg.Id()] = pg;
    return *this;
}

ReadGroupInfo BamHeader::ReadGroup(const std::string& id) const
{
    const auto iter = d_->readGroups_.find(id);
    if (iter == d_->readGroups_.cend()) throw std::runtime_error("read group ID not found");
    return iter->second;
}

std::vector<std::string> BamHeader::ReadGroupIds() const
{
    std::vector<std::string> result;
    result.reserve(d_->readGroups_.size());
    const auto end = d_->readGroups_.cend();
    auto iter = d_->readGroups_.cbegin();
    for (; iter != end; ++iter)
        result.push_back(iter->first);
    return result;
}

std::vector<ReadGroupInfo> BamHeader::ReadGroups() const
{
    std::vector<ReadGroupInfo> result;
    result.reserve(d_->readGroups_.size());
    const auto end = d_->readGroups_.cend();
    auto iter = d_->readGroups_.cbegin();
    for (; iter != end; ++iter)
        result.push_back(iter->second);
    return result;
}

BamHeader& BamHeader::ReadGroups(const std::vector<ReadGroupInfo>& readGroups)
{
    d_->readGroups_.clear();
    for (const auto& rg : readGroups)
        d_->readGroups_[rg.Id()] = rg;
    return *this;
}

SequenceInfo BamHeader::Sequence(const std::string& name) const
{
    // TODO: SequenceId(name) throws if not found, should we do so here as well?

    const auto iter = d_->sequenceIdLookup_.find(name);
    if (iter == d_->sequenceIdLookup_.cend()) return SequenceInfo();
    const auto index = iter->second;
    assert(index >= 0 && (size_t)index < d_->sequences_.size());
    return d_->sequences_.at(index);
}

int32_t BamHeader::SequenceId(const std::string& name) const
{
    const auto iter = d_->sequenceIdLookup_.find(name);
    if (iter == d_->sequenceIdLookup_.cend()) throw std::runtime_error("sequence not found");
    return iter->second;
}

std::vector<std::string> BamHeader::SequenceNames() const
{
    std::vector<std::string> result;
    result.reserve(d_->sequences_.size());
    const auto end = d_->sequences_.cend();
    auto iter = d_->sequences_.cbegin();
    for (; iter != end; ++iter)
        result.push_back(iter->Name());
    return result;
}

BamHeader& BamHeader::Sequences(const std::vector<SequenceInfo>& sequences)
{
    d_->sequences_.clear();
    for (const auto& seq : sequences)
        AddSequence(seq);
    return *this;
}

std::string BamHeader::ToSam() const
{
    // init stream
    std::stringstream out("");

    // @HD
    const auto& outputVersion = (d_->version_.empty() ? std::string(hts_version()) : d_->version_);
    const auto& outputSortOrder =
        (d_->sortOrder_.empty() ? std::string("unknown") : d_->sortOrder_);
    const auto& outputPbBamVersion =
        (d_->pacbioBamVersion_.empty() ? internal::Version::Current.ToString()
                                       : d_->pacbioBamVersion_);

    out << internal::BamHeaderPrefixHD
        << internal::MakeSamTag(internal::BamHeaderTokenVN, outputVersion)
        << internal::MakeSamTag(internal::BamHeaderTokenSO, outputSortOrder)
        << internal::MakeSamTag(internal::BamHeaderTokenpb, outputPbBamVersion) << std::endl;

    // @SQ
    for (const auto& seq : d_->sequences_)
        out << seq.ToSam() << std::endl;

    // @RG
    for (const auto& rgIter : d_->readGroups_)
        out << rgIter.second.ToSam() << std::endl;

    // @PG
    for (const auto& progIter : d_->programs_)
        out << progIter.second.ToSam() << std::endl;

    // @CO
    for (const auto& comment : d_->comments_)
        out << internal::BamHeaderPrefixCO << '\t' << comment << std::endl;

    // return result
    return out.str();
}

}  // namespace BAM
}  // namespace PacBio
