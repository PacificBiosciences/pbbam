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

#include "pbbam/BamHeader.h"
#include "StringUtils.h"
#include "Version.h"
#include <htslib/hts.h>
#include <sstream>
#include <set>
#include <cassert>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

static const string prefix_HD = string("@HD");
static const string prefix_SQ = string("@SQ");
static const string prefix_RG = string("@RG");
static const string prefix_PG = string("@PG");
static const string prefix_CO = string("@CO");

static const string token_VN = string("VN");
static const string token_SO = string("SO");
static const string token_pb = string("pb");

static
void EnsureCanMerge(const BamHeader& lhs, const BamHeader& rhs)
{
    // check compatibility
    const bool samVersionOk = lhs.Version() == rhs.Version();
    const bool sortOrderOk  = lhs.SortOrder() == rhs.SortOrder();
    const bool pbVersionOk  = lhs.PacBioBamVersion() == rhs.PacBioBamVersion();
    const bool sequencesOk  = ( (lhs.SortOrder() == "coordinate") ? lhs.Sequences() == rhs.Sequences()
                                                                  : true);

    // if all checks out, return
    if (samVersionOk && sortOrderOk && pbVersionOk && sequencesOk)
        return;

    // else, format error message & throw
    stringstream e;
    e << "could not merge BAM headers:" << endl;

    if (!samVersionOk) {
        e << "  mismatched SAM versions (@HD:VN) : ("
          << lhs.Version() << ", " << rhs.Version()
          << ")" << endl;
    }

    if (!sortOrderOk) {
        e << "  mismatched sort orders (@HD:SO) : ("
          << lhs.SortOrder() << ", " << rhs.SortOrder()
          << ")" << endl;
    }

    if (!pbVersionOk) {
        e << "  mismatched PacBio BAM versions (@HD:pb) : ("
          << lhs.PacBioBamVersion() << ", " << rhs.PacBioBamVersion()
          << ")" << endl;
    }

    if (!sequencesOk)
        e << "  mismatched sequence lists (@SQ entries)" << endl;

    throw std::runtime_error(e.str());
}

} // namespace internal
} // namespace BAM
} // namespace PacBio

BamHeader::BamHeader(const string& samHeaderText)
    : d_(new internal::BamHeaderPrivate)
{
    istringstream s(samHeaderText);
    string line("");
    string firstToken;
    while (getline(s, line)) {

        // skip if line is not long enough to contain true values
        if (line.length() < 5)
            continue;

        // determine token at beginning of line
        firstToken = line.substr(0,3);

        if (firstToken == internal::prefix_HD) {

            // pop off '@HD\t', then split HD lines into tokens
            const vector<string>& tokens = internal::Split(line.substr(4), '\t');
            for (const string& token : tokens) {
                const string& tokenTag   = token.substr(0,2);
                const string& tokenValue = token.substr(3);

                // set header contents
                if      (tokenTag == internal::token_VN) Version(tokenValue);
                else if (tokenTag == internal::token_SO) SortOrder(tokenValue);
                else if (tokenTag == internal::token_pb) PacBioBamVersion(tokenValue);
            }

            // check for required tags
            if (Version().empty())
                Version(string(hts_version()));
        }

        else if (firstToken == internal::prefix_SQ)
            AddSequence(SequenceInfo::FromSam(line));

        else if (firstToken == internal::prefix_RG)
            AddReadGroup(ReadGroupInfo::FromSam(line));

        else if (firstToken == internal::prefix_PG)
            AddProgram(ProgramInfo::FromSam(line));

        else if (firstToken == internal::prefix_CO)
            AddComment(line.substr(4));
    }
}

BamHeader& BamHeader::operator+=(const BamHeader& other)
{
    internal::EnsureCanMerge(*this, other);

    // merge read groups
    for (const auto& rg : other.ReadGroups()) {
        if (!HasReadGroup(rg.Id()))
            AddReadGroup(rg);
    }

    // merge programs
    for (const auto& pg : other.Programs()) {
        if (!HasProgram(pg.Id()))
            AddProgram(pg);
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

BamHeader& BamHeader::ClearSequences(void)
{
    d_->sequenceIdLookup_.clear();
    d_->sequences_.clear();
    return *this;
}

BamHeader BamHeader::DeepCopy(void) const
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
    const auto fileVersion = internal::Version{ version };
    if (fileVersion < internal::Version::Minimum) {
        auto msg  = string{ "invalid PacBio BAM version number" };
             msg += ( "(" + fileVersion.ToString() + ")");
             msg += string{ "is older than the minimum supported version" };
             msg += ( "(" + internal::Version::Minimum.ToString() + ")");
        throw std::runtime_error(msg);
    }
    return *this;
}

ProgramInfo BamHeader::Program(const std::string& id) const
{
    const auto iter = d_->programs_.find(id);
    if (iter == d_->programs_.cend())
        throw std::runtime_error("program ID not found");
    return iter->second;
}

vector<string> BamHeader::ProgramIds(void) const
{
    vector<string> result;
    result.reserve(d_->programs_.size());
    const auto end = d_->programs_.cend();
    auto iter = d_->programs_.cbegin();
    for ( ; iter != end; ++iter )
        result.push_back(iter->first);
    return result;
}

vector<ProgramInfo> BamHeader::Programs(void) const
{
    vector<ProgramInfo> result;
    result.reserve(d_->programs_.size());
    const auto end = d_->programs_.cend();
    auto iter = d_->programs_.cbegin();
    for ( ; iter != end; ++iter )
        result.push_back(iter->second);
    return result;
}

BamHeader& BamHeader::Programs(const vector<ProgramInfo>& programs)
{
    d_->programs_.clear();
    for (const ProgramInfo& pg : programs)
        d_->programs_[pg.Id()] = pg;
    return *this;
}

ReadGroupInfo BamHeader::ReadGroup(const std::string& id) const
{
    const auto iter = d_->readGroups_.find(id);
    if (iter == d_->readGroups_.cend())
        throw std::runtime_error("read group ID not found");
    return iter->second;
}

vector<string> BamHeader::ReadGroupIds(void) const
{
    vector<string> result;
    result.reserve(d_->readGroups_.size());
    const auto end = d_->readGroups_.cend();
    auto iter = d_->readGroups_.cbegin();
    for ( ; iter != end; ++iter )
        result.push_back(iter->first);
    return result;
}

vector<ReadGroupInfo> BamHeader::ReadGroups(void) const
{
    vector<ReadGroupInfo> result;
    result.reserve(d_->readGroups_.size());
    const auto end = d_->readGroups_.cend();
    auto iter = d_->readGroups_.cbegin();
    for ( ; iter != end; ++iter )
        result.push_back(iter->second);
    return result;
}

BamHeader& BamHeader::ReadGroups(const vector<ReadGroupInfo>& readGroups)
{
    d_->readGroups_.clear();
    for (const ReadGroupInfo& rg : readGroups)
        d_->readGroups_[rg.Id()] = rg;
    return *this;
}

SequenceInfo BamHeader::Sequence(const std::string& name) const
{
    // TODO: SequenceId(name) throws if not found, should we do so here as well?

    const auto iter = d_->sequenceIdLookup_.find(name);
    if (iter == d_->sequenceIdLookup_.cend())
        return SequenceInfo();
    const int index = iter->second;
    assert(index >= 0 && (size_t)index < d_->sequences_.size());
    return d_->sequences_.at(index);
}

int32_t BamHeader::SequenceId(const std::string& name) const
{
    const auto iter = d_->sequenceIdLookup_.find(name);
    if (iter == d_->sequenceIdLookup_.cend())
        throw std::runtime_error("sequence not found");
    return iter->second;
}

vector<string> BamHeader::SequenceNames(void) const
{
    vector<string> result;
    result.reserve(d_->sequences_.size());
    const auto end = d_->sequences_.cend();
    auto iter = d_->sequences_.cbegin();
    for ( ; iter != end; ++iter )
        result.push_back(iter->Name());
    return result;
}

BamHeader& BamHeader::Sequences(const vector<SequenceInfo>& sequences)
{
    d_->sequences_.clear();
    for (const SequenceInfo& seq : sequences)
        AddSequence(seq);
    return *this;
}

string BamHeader::ToSam(void) const
{
    // init stream
    stringstream out("");

    // @HD
    const string& outputVersion   = (d_->version_.empty() ? string(hts_version()) : d_->version_);
    const string& outputSortOrder = (d_->sortOrder_.empty() ? string("unknown") : d_->sortOrder_);
    const string& outputPbBamVersion = (d_->pacbioBamVersion_.empty() ? internal::Version::Current.ToString()
                                                                      : d_->pacbioBamVersion_);

    out << internal::prefix_HD
        << internal::MakeSamTag(internal::token_VN, outputVersion)
        << internal::MakeSamTag(internal::token_SO, outputSortOrder)
        << internal::MakeSamTag(internal::token_pb, outputPbBamVersion)
        << endl;

    // @SQ
    for (const SequenceInfo& seq : d_->sequences_)
        out << seq.ToSam() << endl;

    // @RG
    for (const auto& rgIter : d_->readGroups_)
        out << rgIter.second.ToSam() << endl;

    // @PG
    for (const auto& progIter : d_->programs_)
        out << progIter.second.ToSam() << endl;

    // @CO
    for (const string& comment : d_->comments_)
        out << internal::prefix_CO << '\t' << comment << endl;

    // return result
    return out.str();
}
