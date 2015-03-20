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

#ifndef BAMRECORDBUILDER_H
#define BAMRECORDBUILDER_H

#include "pbbam/BamRecord.h"
#include "pbbam/BamHeader.h"
#include "pbbam/Config.h"
#include <memory>
#include <string>

namespace PacBio {
namespace BAM {

class PBBAM_EXPORT BamImplBuilder
{

};


class PBBAM_EXPORT BamRecordBuilder
{
public:
    /// \name Constructors & Related Methods
    /// \{

    BamRecordBuilder(void);
    explicit BamRecordBuilder(const std::shared_ptr<BamHeader>& header);
    BamRecordBuilder(const BamRecord& prototype);
    BamRecordBuilder(const BamRecordBuilder& other);
    BamRecordBuilder(BamRecordBuilder&& other);
    BamRecordBuilder& operator=(const BamRecordBuilder& other);
    BamRecordBuilder& operator=(BamRecordBuilder&& other);
    ~BamRecordBuilder(void);

    /// \}

public:
    /// \name Record-Building
    /// \{

    /// Builds a BamRecord from current builder attributes
    ///
    /// \returns BamRecord object
    BamRecord Build(void) const;

    /// Replaces an existing BamRecord's data with current builder attributes
    ///
    /// \param[out] record resulting record
    /// \returns true if successful
    bool BuildInPlace(BamRecord& record) const;

    /// Resets builder attributes to default values
    void Reset(void);

    /// Resets builder attributes with existing BamRecord data
    ///
    /// \param[in] prototype
    void Reset(const BamRecord& prototype);

    /// Resets builder attributes with existing BamRecord data
    ///
    /// \param[in] prototype
    void Reset(BamRecord&& prototype);

    /// \}

public:

    /// \name Core Attribute Setup
    /// \{

    /// Sets the record's (BAI) index bin ID.
    ///
    /// \param[in] bin BAI index bin ID.
    /// \returns reference to this builder
    BamRecordBuilder& Bin(const uint32_t bin);

    /// Sets this record's alignment flag, using a raw integer.
    ///
    /// \param[in] flag raw alignment flag
    /// \returns reference to this record
    BamRecordBuilder& Flag(const uint32_t flag);

    /// Sets this record's insert size.
    ///
    /// \param[in] iSize insert size
    /// \returns reference to this record
    BamRecordBuilder& InsertSize(const int32_t iSize);

    /// Sets this record's map quality.
    ///
    /// \param[in] mapQual mapping quality - value of 255 indicates "unknown"
    /// \returns reference to this record
    BamRecordBuilder& MapQuality(const uint8_t mapQual);

    /// Sets this record's mate's mapped position.
    ///
    /// \param[in] pos mapped position. A value of -1 indicates unmapped.
    /// \returns reference to this record
    BamRecordBuilder& MatePosition(const int32_t pos);

    /// Sets this record's mate's mapped reference ID
    ///
    /// \param[in] id reference ID. A value of -1 indicates unmapped.
    /// \returns reference to this record
    BamRecordBuilder& MateReferenceId(const int32_t id);

    /// Sets this record's mapped position.
    ///
    /// \param[in] pos mapped position. A value of -1 indicates unmapped.
    /// \returns reference to this record
    BamRecordBuilder& Position(const int32_t pos);

    /// Sets this record's mapped reference ID
    ///
    /// \param[in] id reference ID. A value of -1 indicates unmapped.
    /// \returns reference to this record
    BamRecordBuilder& ReferenceId(const int32_t id);

    /// \}

public:
    /// \name Alignment Flag Setup
    /// \{

    /// Sets whether this record is a PCR/optical duplicate
    BamRecordBuilder& SetDuplicate(bool ok);

    /// Sets whether this record failed quality controls
    BamRecordBuilder& SetFailedQC(bool ok);

    /// Sets whether this record is the first mate of a pair.
    BamRecordBuilder& SetFirstMate(bool ok);

    /// Sets whether this record was aligned.
    BamRecordBuilder& SetMapped(bool ok);

    /// Sets whether this record's mate was aligned.
    BamRecordBuilder& SetMateMapped(bool ok);

    /// Sets whether this record's mate mapped to reverse strand.
    BamRecordBuilder& SetMateReverseStrand(bool ok);

    /// Sets whether this record came from paired-end sequencing.
    BamRecordBuilder& SetPaired(bool ok);

    /// Sets whether this record is a read's primary alignment.
    BamRecordBuilder& SetPrimaryAlignment(bool ok);

    /// Sets whether this record & its mate were properly mapped, per the aligner.
    BamRecordBuilder& SetProperPair(bool ok);

    /// Sets whether this record mapped to reverse strand.
    BamRecordBuilder& SetReverseStrand(bool ok);

    /// Sets whether this record is the second mate of a pair.
    BamRecordBuilder& SetSecondMate(bool ok);

    /// Sets whether this record is a supplementary alignment.
    BamRecordBuilder& SetSupplementaryAlignment(bool ok);

    /// \}

public:
    /// \name Variable-Length Data Setup
    /// \{

    BamRecordBuilder& Name(const std::string& name);
    BamRecordBuilder& Name(std::string&& name);

    BamRecordBuilder& Sequence(const std::string& sequence);
    BamRecordBuilder& Sequence(std::string&& sequence);

    BamRecordBuilder& Qualities(const std::string& qualities);
    BamRecordBuilder& Qualities(std::string&& qualities);

    BamRecordBuilder& Cigar(const PacBio::BAM::Cigar& cigar);
    BamRecordBuilder& Cigar(PacBio::BAM::Cigar&& cigar);

    BamRecordBuilder& Tags(const TagCollection& tags);
    BamRecordBuilder& Tags(TagCollection&& tags);

private:
    std::shared_ptr<BamHeader> header_;

    bam1_core_t core_;
    std::string name_;
    std::string sequence_;
    std::string qualities_;
    PacBio::BAM::Cigar cigar_;
    TagCollection tags_;
};

inline BamRecordBuilder& BamRecordBuilder::Bin(const uint32_t bin)
{ core_.bin = bin; return *this; }

inline BamRecordBuilder& BamRecordBuilder::Flag(const uint32_t flag)
{ core_.flag = flag; return *this; }

inline BamRecordBuilder& BamRecordBuilder::InsertSize(const int32_t iSize)
{ core_.isize = iSize; return *this; }

inline BamRecordBuilder& BamRecordBuilder::MapQuality(const uint8_t mapQual)
{ core_.qual = mapQual; return *this; }

inline BamRecordBuilder& BamRecordBuilder::MatePosition(const int32_t pos)
{ core_.mpos = pos; return *this; }

inline BamRecordBuilder& BamRecordBuilder::MateReferenceId(const int32_t id)
{ core_.mtid = id; return *this; }

inline BamRecordBuilder& BamRecordBuilder::Position(const int32_t pos)
{ core_.pos = pos; return *this; }

inline BamRecordBuilder& BamRecordBuilder::Qualities(const std::string& qualities)
{ qualities_ = qualities; return *this; }

inline BamRecordBuilder& BamRecordBuilder::Qualities(std::string&& qualities)
{ qualities_ = std::move(qualities); return *this; }

inline BamRecordBuilder& BamRecordBuilder::ReferenceId(const int32_t id)
{ core_.tid = id; return *this; }

inline BamRecordBuilder& BamRecordBuilder::Tags(const TagCollection& tags)
{ tags_ = tags; return *this; }

inline BamRecordBuilder& BamRecordBuilder::Tags(TagCollection&& tags)
{ tags_ = std::move(tags); return *this; }

} // namespace BAM
} // namespace PacBio

#endif // BAMRECORDBUILDER_H
