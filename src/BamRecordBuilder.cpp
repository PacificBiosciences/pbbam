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

#include "pbbam/BamRecordBuilder.h"
#include "pbbam/BamTagCodec.h"
#include "AssertUtils.h"
#include "MemoryUtils.h"
#include <htslib/sam.h>
#include <cstring>
#include <memory>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

BamRecordBuilder::BamRecordBuilder(void)
{
    // ensure proper clean slate
    Reset();

    // initialize with some space for data
    name_.reserve(256);
    sequence_.reserve(2096);
    qualities_.reserve(2096);
    cigar_.reserve(256);
}

BamRecordBuilder::BamRecordBuilder(const BamHeader& header)
    : header_(header)
{
    // ensure proper clean slate
    Reset();

    // initialize with some space for data
    name_.reserve(256);
    sequence_.reserve(2096);
    qualities_.reserve(2096);
    cigar_.reserve(256);
}

BamRecordBuilder::BamRecordBuilder(const BamRecord& prototype)
    : header_(prototype.Header())
{
    Reset(prototype);
}

BamRecordBuilder::BamRecordBuilder(const BamRecordBuilder& other)
    : core_(other.core_)
    , name_(other.name_)
    , sequence_(other.sequence_)
    , qualities_(other.qualities_)
    , cigar_(other.cigar_)
    , tags_(other.tags_)
{ }

BamRecordBuilder::BamRecordBuilder(BamRecordBuilder&& other)
    : core_(std::move(other.core_))
    , name_(std::move(other.name_))
    , sequence_(std::move(other.sequence_))
    , qualities_(std::move(other.qualities_))
    , cigar_(std::move(other.cigar_))
    , tags_(std::move(other.tags_))
{  }

BamRecordBuilder& BamRecordBuilder::operator=(const BamRecordBuilder& other)
{
    core_ = other.core_;
    name_ = other.name_;
    sequence_  = other.sequence_;
    qualities_ = other.qualities_;
    cigar_ = other.cigar_;
    tags_  = other.tags_;
    return *this;
}

BamRecordBuilder& BamRecordBuilder::operator=(BamRecordBuilder&& other)
{
    core_ = std::move(other.core_);
    name_ = std::move(other.name_);
    sequence_  = std::move(other.sequence_);
    qualities_ = std::move(other.qualities_);
    cigar_ = std::move(other.cigar_);
    tags_  = std::move(other.tags_);
    return *this;
}

BamRecordBuilder::~BamRecordBuilder(void) { }

BamRecord BamRecordBuilder::Build(void) const
{
    BamRecord result(header_);
    BuildInPlace(result);
    return result;
}

bool BamRecordBuilder::BuildInPlace(BamRecord& record) const
{
    // initialize with basic 'core data'
    PBBAM_SHARED_PTR<bam1_t> recordRawData = internal::BamRecordMemory::GetRawData(record); /*   record.impl_.RawData().get();*/
    PB_ASSERT_OR_RETURN_VALUE(recordRawData, false);
    PB_ASSERT_OR_RETURN_VALUE(recordRawData->data, false);
    recordRawData->core = core_;

    // setup variable length data
    const vector<uint8_t> encodedTags = std::move(BamTagCodec::Encode(tags_));

    const size_t nameLength  = name_.size() + 1;
    const size_t numCigarOps = cigar_.size();
    const size_t cigarLength = numCigarOps * sizeof(uint32_t);
    const size_t seqLength   = sequence_.size();
    const size_t qualLength  = seqLength;
    const size_t tagLength   = encodedTags.size();
    const size_t dataLength  = nameLength + cigarLength + seqLength + qualLength + tagLength;

    // realloc if necessary
    uint8_t* varLengthDataBlock = recordRawData->data;
    PB_ASSERT_OR_RETURN_VALUE(varLengthDataBlock, false);
    size_t allocatedDataLength = recordRawData->m_data;
    if (allocatedDataLength < dataLength) {
        allocatedDataLength = dataLength;
        kroundup32(allocatedDataLength);
        varLengthDataBlock = (uint8_t*)realloc(varLengthDataBlock, allocatedDataLength);
    }
    recordRawData->data = varLengthDataBlock;
    recordRawData->l_data = dataLength;
    recordRawData->m_data = allocatedDataLength;

    size_t index = 0;

    // name
    memcpy(&varLengthDataBlock[index], name_.c_str(), nameLength);
    index += nameLength;

    // cigar
    if (cigarLength > 0) {
        vector<uint32_t> encodedCigar(numCigarOps);
        for (size_t i = 0; i < numCigarOps; ++i) {
            const CigarOperation& op = cigar_.at(i);
            encodedCigar[i] = op.Length() << BAM_CIGAR_SHIFT;
            const uint8_t type = static_cast<uint8_t>(op.Type());
            PB_ASSERT_OR_RETURN_VALUE(type >= 0 && type < 8, false);
            encodedCigar[i] |= type;
        }
        memcpy(&varLengthDataBlock[index], &encodedCigar[0], cigarLength);
        index += cigarLength;

        // update bin after we've calculated cigar info
        const int32_t endPosition = bam_cigar2rlen(recordRawData->core.n_cigar, &encodedCigar[0]);
        recordRawData->core.bin = hts_reg2bin(core_.pos, endPosition, 14, 5);
    }

    // seq & qual
    if (seqLength > 0) {

        uint8_t* s = &varLengthDataBlock[index];
        for (size_t i = 0; i < seqLength; ++i)
            s[i>>1] |= ( seq_nt16_table[static_cast<int>(sequence_.at(i))] << ((~i&1)<<2) );
        index += seqLength;

        uint8_t* q = &varLengthDataBlock[index];
        if (!qualities_.empty())
            memset(q, 0xFF, seqLength);
        else {
            for (size_t i = 0; i < seqLength; ++i)
                q[i] = qualities_.at(i) - 33;
        }
        index += seqLength;
    }

    // tags
    if (tagLength > 0) {
        PB_ASSERT_OR_RETURN_VALUE(!encodedTags.empty(), false);
        memcpy(&varLengthDataBlock[index], &encodedTags[0], tagLength);
        index += tagLength;
    }

    // sanity check
    PB_ASSERT_OR_RETURN_VALUE(index == dataLength, false);
    return true;
}

BamRecordBuilder& BamRecordBuilder::Cigar(const PacBio::BAM::Cigar& cigar)
{
    core_.n_cigar = cigar.size();
    cigar_ = cigar;
    return *this;
}

BamRecordBuilder& BamRecordBuilder::Cigar(PacBio::BAM::Cigar&& cigar)
{
    core_.n_cigar = cigar.size();
    cigar_ = std::move(cigar);
    return *this;
}

BamRecordBuilder& BamRecordBuilder::Name(const std::string& name)
{
    core_.l_qname = name.size() + 1; // (NULL-term)
    name_ = name;
    return *this;
}

BamRecordBuilder& BamRecordBuilder::Name(std::string&& name)
{
    core_.l_qname = name.size() + 1; // (NULL-term)
    name_ = std::move(name);
    return *this;
}

void BamRecordBuilder::Reset(void)
{
    // zeroize fixed-length data
    memset(&core_, 0, sizeof(bam1_core_t));
    core_.l_qname = 1; // always has a NULL-term

    // reset variable-length data
    name_.clear();
    sequence_.clear();
    qualities_.clear();
    cigar_.clear();
    tags_.clear();
}

void BamRecordBuilder::Reset(const BamRecord& prototype)
{
    // ensure clean slate
    Reset();
    header_ = prototype.Header();

    // reset core data
    const PBBAM_SHARED_PTR<bam1_t> rawData = internal::BamRecordMemory::GetRawData(prototype); //  prototype.impl_.RawData().get();
    PB_ASSERT_OR_RETURN(rawData);
    core_ = rawData->core;

    // reset variable-length data
    const BamRecordImpl& impl = internal::BamRecordMemory::GetImpl(prototype);
    name_ = impl.Name();
    sequence_ = impl.Sequence();
    qualities_ = impl.Qualities().Fastq();
    cigar_ = impl.CigarData();
    tags_ = impl.Tags();
}

void BamRecordBuilder::Reset(BamRecord&& prototype)
{
    // ensure clean slate
    Reset();
    header_ = std::move(prototype.Header());

    // reset core data
    const PBBAM_SHARED_PTR<bam1_t> rawData = internal::BamRecordMemory::GetRawData(prototype); //  prototype.impl_.RawData().get();
    PB_ASSERT_OR_RETURN(rawData);
    core_ = std::move(rawData->core);

    // reset variable-length data
    const BamRecordImpl& impl = internal::BamRecordMemory::GetImpl(prototype);
    name_ = std::move(impl.Name());
    sequence_ = std::move(impl.Sequence());
    qualities_ = std::move(impl.Qualities().Fastq());
    cigar_ = std::move(impl.CigarData());
    tags_ = std::move(impl.Tags());
}

BamRecordBuilder& BamRecordBuilder::Sequence(const std::string& sequence)
{
    core_.l_qseq = sequence.size();
    sequence_ = sequence;
    return *this;
}

BamRecordBuilder& BamRecordBuilder::Sequence(std::string&& sequence)
{
    core_.l_qseq = sequence.size();
    sequence_ = std::move(sequence);
    return *this;
}

BamRecordBuilder& BamRecordBuilder::SetDuplicate(bool ok)
{
    if (ok) core_.flag |=  BamRecordImpl::DUPLICATE;
    else    core_.flag &= ~BamRecordImpl::DUPLICATE;
    return *this;
}

BamRecordBuilder& BamRecordBuilder::SetFailedQC(bool ok)
{
    if (ok) core_.flag |=  BamRecordImpl::FAILED_QC;
    else    core_.flag &= ~BamRecordImpl::FAILED_QC;
    return *this;
}

BamRecordBuilder& BamRecordBuilder::SetFirstMate(bool ok)
{
    if (ok) core_.flag |=  BamRecordImpl::MATE_1;
    else    core_.flag &= ~BamRecordImpl::MATE_1;
    return *this;
}

BamRecordBuilder& BamRecordBuilder::SetMapped(bool ok)
{
    if (ok) core_.flag &= ~BamRecordImpl::UNMAPPED;
    else    core_.flag |=  BamRecordImpl::UNMAPPED;
    return *this;
}

BamRecordBuilder& BamRecordBuilder::SetMateMapped(bool ok)
{
    if (ok) core_.flag &= ~BamRecordImpl::MATE_UNMAPPED;
    else    core_.flag |=  BamRecordImpl::MATE_UNMAPPED;
    return *this;
}

BamRecordBuilder& BamRecordBuilder::SetMateReverseStrand(bool ok)
{
    if (ok) core_.flag |=  BamRecordImpl::MATE_REVERSE_STRAND;
    else    core_.flag &= ~BamRecordImpl::MATE_REVERSE_STRAND;
    return *this;
}

BamRecordBuilder& BamRecordBuilder::SetPaired(bool ok)
{
    if (ok) core_.flag |=  BamRecordImpl::PAIRED;
    else    core_.flag &= ~BamRecordImpl::PAIRED;
    return *this;
}

BamRecordBuilder& BamRecordBuilder::SetPrimaryAlignment(bool ok)
{
    if (ok) core_.flag &= ~BamRecordImpl::SECONDARY;
    else    core_.flag |=  BamRecordImpl::SECONDARY;
    return *this;
}

BamRecordBuilder& BamRecordBuilder::SetProperPair(bool ok)
{
    if (ok) core_.flag |=  BamRecordImpl::PROPER_PAIR;
    else    core_.flag &= ~BamRecordImpl::PROPER_PAIR;
    return *this;
}

BamRecordBuilder& BamRecordBuilder::SetReverseStrand(bool ok)
{
    if (ok) core_.flag |=  BamRecordImpl::REVERSE_STRAND;
    else    core_.flag &= ~BamRecordImpl::REVERSE_STRAND;
    return *this;
}

BamRecordBuilder& BamRecordBuilder::SetSecondMate(bool ok)
{
    if (ok) core_.flag |=  BamRecordImpl::MATE_2;
    else    core_.flag &= ~BamRecordImpl::MATE_2;
    return *this;
}

BamRecordBuilder& BamRecordBuilder::SetSupplementaryAlignment(bool ok)
{
    if (ok) core_.flag |=  BamRecordImpl::SUPPLEMENTARY;
    else    core_.flag &= ~BamRecordImpl::SUPPLEMENTARY;
    return *this;
}
