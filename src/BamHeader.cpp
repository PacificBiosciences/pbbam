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

#include "pbbam/BamHeader.h"
#include "SequenceUtils.h"
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

} // namespace internal
} // namespace BAM
} // namespace PacBio

BamHeader::BamHeader(void) { }

BamHeader::BamHeader(const BamHeader& other)
    : version_(other.version_)
    , pacbioBamVersion_(other.pacbioBamVersion_)
    , sortOrder_(other.sortOrder_)
    , readGroups_(other.readGroups_)
    , programs_(other.programs_)
    , comments_(other.comments_)
    , sequences_(other.sequences_)
    , sequenceIdLookup_(other.sequenceIdLookup_)
{ }

BamHeader::BamHeader(BamHeader&& other)
    : version_(std::move(other.version_))
    , pacbioBamVersion_(std::move(other.pacbioBamVersion_))
    , sortOrder_(std::move(other.sortOrder_))
    , readGroups_(std::move(other.readGroups_))
    , programs_(std::move(other.programs_))
    , comments_(std::move(other.comments_))
    , sequences_(std::move(other.sequences_))
    , sequenceIdLookup_(std::move(other.sequenceIdLookup_))
{ }

BamHeader::~BamHeader(void) { }

BamHeader& BamHeader::operator=(const BamHeader& other)
{
    version_ = other.version_;
    pacbioBamVersion_ = other.pacbioBamVersion_;
    sortOrder_ = other.sortOrder_;
    readGroups_ = other.readGroups_;
    programs_ = other.programs_;
    comments_ = other.comments_;
    sequences_ = other.sequences_;
    sequenceIdLookup_ = other.sequenceIdLookup_;
    return *this;
}

BamHeader& BamHeader::operator=(BamHeader&& other)
{
    version_ = std::move(other.version_);
    pacbioBamVersion_ = std::move(other.pacbioBamVersion_);
    sortOrder_ = std::move(other.sortOrder_);
    readGroups_ = std::move(other.readGroups_);
    programs_ = std::move(other.programs_);
    comments_ = std::move(other.comments_);
    sequences_ = std::move(other.sequences_);
    sequenceIdLookup_ = std::move(other.sequenceIdLookup_);
    return *this;
}

BamHeader& BamHeader::AddSequence(const SequenceInfo& sequence)
{
    sequences_.push_back(sequence);
    sequenceIdLookup_[sequence.Name()] = sequences_.size() - 1;
    return *this;
}

std::shared_ptr<BamHeader> BamHeader::FromSam(const string& sam)
{
    std::shared_ptr<BamHeader> header(new BamHeader);

    istringstream s(sam);
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
                if      (tokenTag == internal::token_VN) header->Version(tokenValue);
                else if (tokenTag == internal::token_SO) header->SortOrder(tokenValue);
                else if (tokenTag == internal::token_pb) header->PacBioBamVersion(tokenValue);
            }

            // check for required tags
            if (header->Version().empty())
                header->Version(string(hts_version()));
        }

        else if (firstToken == internal::prefix_SQ)
            header->AddSequence(SequenceInfo::FromSam(line));

        else if (firstToken == internal::prefix_RG)
            header->AddReadGroup(ReadGroupInfo::FromSam(line));

        else if (firstToken == internal::prefix_PG)
            header->AddProgram(ProgramInfo::FromSam(line));

        else if (firstToken == internal::prefix_CO)
            header->AddComment(line.substr(4));
    }

    return header;
}

vector<string> BamHeader::ProgramIds(void) const
{
    vector<string> result;
    result.reserve(programs_.size());
    const auto end = programs_.cend();
    auto iter = programs_.cbegin();
    for ( ; iter != end; ++iter )
        result.push_back(iter->first);
    return result;
}

vector<ProgramInfo> BamHeader::Programs(void) const
{
    vector<ProgramInfo> result;
    result.reserve(programs_.size());
    const auto end = programs_.cend();
    auto iter = programs_.cbegin();
    for ( ; iter != end; ++iter )
        result.push_back(iter->second);
    return result;
}

BamHeader& BamHeader::Programs(const vector<ProgramInfo>& programs)
{
    programs_.clear();
    for (const ProgramInfo& pg : programs)
        programs_[pg.Id()] = pg;
    return *this;
}

vector<string> BamHeader::ReadGroupIds(void) const
{
    vector<string> result;
    result.reserve(readGroups_.size());
    const auto end = readGroups_.cend();
    auto iter = readGroups_.cbegin();
    for ( ; iter != end; ++iter )
        result.push_back(iter->first);
    return result;
}

vector<ReadGroupInfo> BamHeader::ReadGroups(void) const
{
    vector<ReadGroupInfo> result;
    result.reserve(readGroups_.size());
    const auto end = readGroups_.cend();
    auto iter = readGroups_.cbegin();
    for ( ; iter != end; ++iter )
        result.push_back(iter->second);
    return result;
}

BamHeader& BamHeader::ReadGroups(const vector<ReadGroupInfo>& readGroups)
{
    readGroups_.clear();
    for (const ReadGroupInfo& rg : readGroups)
        readGroups_[rg.Id()] = rg;
    return *this;
}

SequenceInfo BamHeader::Sequence(const std::string& name) const
{
    const auto iter = sequenceIdLookup_.find(name);
    if (iter == sequenceIdLookup_.cend())
        return SequenceInfo();
    const int index = iter->second;
    assert(index >= 0 && (size_t)index < sequences_.size());
    return sequences_.at(index);
}

vector<string> BamHeader::SequenceNames(void) const
{
    vector<string> result;
    result.reserve(sequences_.size());
    const auto end = sequences_.cend();
    auto iter = sequences_.cbegin();
    for ( ; iter != end; ++iter )
        result.push_back(iter->Name());
    return result;
}

BamHeader& BamHeader::Sequences(const vector<SequenceInfo>& sequences)
{
    sequences_.clear();
    for (const SequenceInfo& seq : sequences)
        AddSequence(seq);
    return *this;
}

string BamHeader::ToSam(void) const
{
    // clear out stream
    stringstream out("");

    // @HD
    const string& outputVersion   = (version_.empty()   ? string(hts_version()) : version_);
    const string& outputSortOrder = (sortOrder_.empty() ? string("unknown") : sortOrder_);
//    const string& outputPbBamVersion = (pbVersion.empty() ? string("3.0b3") : pbVersion);

    out << internal::prefix_HD
        << internal::MakeSamTag(internal::token_VN, outputVersion)
        << internal::MakeSamTag(internal::token_SO, outputSortOrder);
    if (!pacbioBamVersion_.empty())
        out << internal::MakeSamTag(internal::token_pb, pacbioBamVersion_);
     out << endl;

    // @SQ
    for (const SequenceInfo& seq : sequences_)
        out << seq.ToSam() << endl;

    // @RG
    for (const auto& rgIter : readGroups_)
        out << rgIter.second.ToSam() << endl;

    // @PG
    for (const auto& progIter : programs_)
        out  << progIter.second.ToSam() << endl;

    // @CO
    for (const string& comment : comments_)
        out << internal::prefix_CO << '\t' << comment << endl;

    // return result
    return out.str();
}
