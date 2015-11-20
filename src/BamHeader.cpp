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

struct PacBioVersion
{
public:
    PacBioVersion(int major, int minor, int revision)
        : major_(major)
        , minor_(minor)
        , revision_(revision)
    { }

    PacBioVersion(const string& v)
        : major_(0)
        , minor_(0)
        , revision_(0)
    {
        if (v.empty()) {
            auto msg = string{ "PacBio BAM version number is missing (@HD pb:<version> tag). See spec for details." };
            throw std::runtime_error(msg);
        }

        if (v.find('b') != string::npos) {
            auto msg = string{ "invalid version number (" + v + "): beta version BAMs are no longer supported" };
            throw std::runtime_error(msg);
        }

        try {
            const auto fields = Split(v, '.');
            const auto numFields = fields.size();
            if (numFields > 0) {
                major_ = stoi(fields.at(0));
                if (numFields > 1) {
                    minor_ = stoi(fields.at(1));
                    if (numFields > 2 )
                        revision_ = stoi(fields.at(2));
                }
            }
        } catch (std::exception&) {
            auto msg = string{ "invalid version number (" + v + "): failed to parse" };
            throw std::runtime_error(msg);
        }
    }

public:
    bool operator==(const PacBioVersion& other) const
    {
        return major_ == other.major_ &&
               minor_ == other.minor_ &&
               revision_ == other.revision_;
    }

    bool operator<(const PacBioVersion& other) const
    {
        // 2.* < 3.*
        if (major_ < other.major_)
            return true;

        // 3. ==  3.
        else if (major_ == other.major_) {

            // 3.1.* < 3.2.*
            if (minor_ < other.minor_)
                return true;

            // 3.2. == 3.2.
            else if (minor_ == other.minor_) {

                // 3.2.1 < 3.2.2
                if (revision_ < other.revision_)
                    return true;
            }
        }

        // otherwise not less-than
        return false;
    }
    bool operator>=(const PacBioVersion& other) const
    { return !operator<(other); }

public:
    string ToString(void) const
    {
        stringstream s;
        s << major_ << '.' << minor_ << '.' << revision_;
        return s.str();
    }

    string ToMsgString(void) const
    {
        stringstream s;
        s << '(' << ToString() << ')';
        return s.str();
    }

private:
    int major_;
    int minor_;
    int revision_;
};

static const PacBioVersion minimum_version = PacBioVersion(3,0,1);
static const PacBioVersion current_version = PacBioVersion(3,0,1);

} // namespace internal
} // namespace BAM
} // namespace PacBio

BamHeader::BamHeader(void)
    : d_(new internal::BamHeaderPrivate)
{ }

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

BamHeader::BamHeader(const BamHeader& other)
    : d_(other.d_)
{ }

BamHeader::BamHeader(BamHeader&& other)
    : d_(std::move(other.d_))
{ }

BamHeader& BamHeader::operator=(const BamHeader& other)
{ d_ = other.d_; return *this; }

BamHeader& BamHeader::operator=(BamHeader&& other)
{ d_ = std::move(other.d_); return *this; }

BamHeader::~BamHeader(void) { }

BamHeader& BamHeader::AddComment(const std::string& comment)
{ d_->comments_.push_back(comment); return *this; }

BamHeader& BamHeader::AddProgram(const ProgramInfo& pg)
{ d_->programs_[pg.Id()] = pg; return *this; }

BamHeader& BamHeader::AddReadGroup(const ReadGroupInfo& readGroup)
{ d_->readGroups_[readGroup.Id()] = readGroup; return *this; }

BamHeader& BamHeader::AddSequence(const SequenceInfo& sequence)
{
    d_->sequences_.push_back(sequence);
    d_->sequenceIdLookup_[sequence.Name()] = d_->sequences_.size() - 1;
    return *this;
}

BamHeader& BamHeader::ClearComments(void)
{ d_->comments_.clear(); return* this; }

BamHeader& BamHeader::ClearPrograms(void)
{ d_->programs_.clear(); return *this; }

BamHeader& BamHeader::ClearReadGroups(void)
{ d_->readGroups_.clear(); return *this; }

BamHeader& BamHeader::ClearSequences(void)
{
    d_->sequenceIdLookup_.clear();
    d_->sequences_.clear();
    return *this;
}

std::vector<std::string> BamHeader::Comments(void) const
{ return d_->comments_; }

BamHeader& BamHeader::Comments(const std::vector<std::string>& comments)
{ d_->comments_ = comments; return *this; }

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

bool BamHeader::HasProgram(const std::string& id) const
{ return d_->programs_.find(id) != d_->programs_.cend(); }

bool BamHeader::HasReadGroup(const std::string& id) const
{ return d_->readGroups_.find(id) != d_->readGroups_.cend(); }

bool BamHeader::HasSequence(const std::string& name) const
{ return d_->sequenceIdLookup_.find(name) != d_->sequenceIdLookup_.cend(); }

std::string BamHeader::PacBioBamVersion(void) const
{ return d_->pacbioBamVersion_; }

BamHeader& BamHeader::PacBioBamVersion(const std::string& version)
{
    const auto fileVersion = internal::PacBioVersion{ version };
    if (fileVersion >= internal::minimum_version)
        d_->pacbioBamVersion_ = version;
    else {
        d_->pacbioBamVersion_.clear();
        auto msg  = string{ "invalid PacBio BAM version number" };
             msg += fileVersion.ToMsgString();
             msg += string{ "is older than the minimum supported version" };
             msg += internal::minimum_version.ToMsgString();
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

SequenceInfo BamHeader::Sequence(const int32_t id) const
{
    // throws out of range
    return d_->sequences_.at(id);
}

size_t BamHeader::NumSequences(void) const
{ return d_->sequences_.size(); }

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

std::string BamHeader::SequenceLength(const int32_t id) const
{ return Sequence(id).Length(); }

std::string BamHeader::SequenceName(const int32_t id) const
{ return Sequence(id).Name(); }

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

std::vector<SequenceInfo> BamHeader::Sequences(void) const
{ return d_->sequences_; }

BamHeader& BamHeader::Sequences(const vector<SequenceInfo>& sequences)
{
    d_->sequences_.clear();
    for (const SequenceInfo& seq : sequences)
        AddSequence(seq);
    return *this;
}

std::string BamHeader::SortOrder(void) const
{ return d_->sortOrder_; }

BamHeader& BamHeader::SortOrder(const std::string& order)
{ d_->sortOrder_ = order; return *this; }

string BamHeader::ToSam(void) const
{
    // clear out stream
    stringstream out("");

    // @HD
    const string& outputVersion   = (d_->version_.empty()   ? string(hts_version()) : d_->version_);
    const string& outputSortOrder = (d_->sortOrder_.empty() ? string("unknown") : d_->sortOrder_);
    const string& outputPbBamVersion = (d_->pacbioBamVersion_.empty() ? internal::current_version.ToString()
                                                                      : d_->pacbioBamVersion_);

    out << internal::prefix_HD
        << internal::MakeSamTag(internal::token_VN, outputVersion)
        << internal::MakeSamTag(internal::token_SO, outputSortOrder)
        << internal::MakeSamTag(internal::token_pb, outputPbBamVersion)
        << endl;

//    if (!d_->pacbioBamVersion_.empty())
//        out << internal::MakeSamTag(internal::token_pb, d_->pacbioBamVersion_);
//     out << endl;

    // @SQ
    for (const SequenceInfo& seq : d_->sequences_)
        out << seq.ToSam() << endl;

    // @RG
    for (const auto& rgIter : d_->readGroups_)
        out << rgIter.second.ToSam() << endl;

    // @PG
    for (const auto& progIter : d_->programs_)
        out  << progIter.second.ToSam() << endl;

    // @CO
    for (const string& comment : d_->comments_)
        out << internal::prefix_CO << '\t' << comment << endl;

    // return result
    return out.str();
}

std::string BamHeader::Version(void) const
{ return d_->version_; }

BamHeader& BamHeader::Version(const std::string& version)
{ d_->version_ = version; return *this; }
