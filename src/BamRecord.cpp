// Copyright (c) 2014, Pacific Biosciences of California, Inc.
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

#include "pbbam/BamRecord.h"
#include "pbbam/BamTagCodec.h"
#include "AssertUtils.h"
#include "MemoryUtils.h"
#include <algorithm>
#include <iostream>
#include <utility>
#include <cstdlib>
#include <cstring>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

BamRecord::BamRecord(void)
    : d_(nullptr)
{
    InitializeData();
}

BamRecord::BamRecord(const BamRecord& other)
    : d_(bam_dup1(other.d_.get()), internal::HtslibRecordDeleter())
{ }

BamRecord::BamRecord(BamRecord&& other)
    : d_(nullptr)
{
    d_.swap(other.d_);
    other.d_.reset();
}

BamRecord& BamRecord::operator=(const BamRecord& other)
{
    if (this != & other) {
        if (d_ == nullptr)
            InitializeData();
        bam_copy1(d_.get(), other.d_.get());
    }
    return *this;
}

BamRecord& BamRecord::operator=(BamRecord&& other)
{
    if (this != & other) {
        d_.swap(other.d_);
        other.d_.reset();
    }
    return *this;
}

BamRecord::~BamRecord(void) { }

bool BamRecord::AddTag(const string& tagName, const Tag& value)
{
    if (tagName.size() != 2 || HasTag(tagName))
        return false;

    const vector<uint8_t>& rawData = BamTagCodec::ToRawData(value);
    if (rawData.empty())
        return false;

    bam_aux_append(d_.get(),
                   tagName.c_str(),
                   BamTagCodec::TagTypeCode(value),
                   rawData.size(),
                   const_cast<uint8_t*>(rawData.data()));
    return true;
}

Cigar BamRecord::CigarData(void) const
{
    Cigar result;
    result.reserve(d_->core.n_cigar);
    uint32_t* cigarData = bam_get_cigar(d_);
    for (uint32_t i = 0; i < d_->core.n_cigar; ++i) {
        const uint32_t length = bam_cigar_oplen(cigarData[i]);
        const char type = bam_cigar_opchr(cigarData[i]);
        result.push_back(CigarOperation(type, length));
    }

    return result;
}

BamRecord& BamRecord::CigarData(const Cigar& cigar)
{
    // determine change in memory needed
    // diffNumBytes: pos -> growing, neg -> shrinking
    const size_t numCigarOps = cigar.size();
    const int diffNumCigars = numCigarOps - d_->core.n_cigar;
    const int diffNumBytes  = diffNumCigars * sizeof(uint32_t);
    const int oldLengthData = d_->l_data;
    d_->l_data += diffNumBytes;
    MaybeReallocData();

    // shift trailing data (seq, qual, tags) as needed
    const uint8_t* oldSequenceStart = bam_get_seq(d_);
    const size_t trailingDataLength = oldLengthData - (oldSequenceStart - d_->data);
    d_->core.n_cigar = numCigarOps;
    uint8_t* newSequenceStart = bam_get_seq(d_);
    memmove(newSequenceStart, oldSequenceStart, trailingDataLength);

    // fill in new CIGAR data
    uint32_t* cigarDataStart = bam_get_cigar(d_);
    for (size_t i = 0; i < numCigarOps; ++i) {
        const CigarOperation& cigarOp = cigar.at(i);
        cigarDataStart[i] = bam_cigar_gen(cigarOp.Length(), static_cast<int>(cigarOp.Operation()));
    }

    return *this;
}

BamRecord& BamRecord::CigarData(const std::string& cigarString)
{
    return CigarData(Cigar::FromStdString(cigarString));
}

bool BamRecord::EditTag(const string& tagName, const Tag& newValue)
{
    return RemoveTag(tagName) && AddTag(tagName, newValue);
}

BamRecord BamRecord::FromRawData(const shared_ptr<bam1_t>& rawData)
{
    BamRecord result;
    bam_copy1(result.d_.get(), rawData.get());
    return result;
}

bool BamRecord::HasTag(const string& tagName) const
{
    if (tagName.size() != 2)
        return false;
    return bam_aux_get(d_.get(), tagName.c_str()) != 0;
}

void BamRecord::InitializeData(void)
{
    d_.reset(bam_init1(), internal::HtslibRecordDeleter());
    d_->data = static_cast<uint8_t*>(calloc(0x800, sizeof(uint8_t)));   // maybe make this value tune-able later?

    // initialized with NULL term for qname
    d_->core.l_qname = 1;
    d_->l_data = 1;
    d_->m_data = 0x800;
}

void BamRecord::MaybeReallocData(void)
{
    // about to grow data contents to l_data size, but m_data is our current max.
    // so we may need to grow. if so, use kroundup to double to next power of 2
    if (d_->m_data < d_->l_data) {
        d_->m_data = d_->l_data;
        kroundup32(d_->m_data);
        d_->data = static_cast<uint8_t*>(realloc(d_->data, d_->m_data));
    }
}

string BamRecord::Name(void) const
{
    return string(bam_get_qname(d_));
}

BamRecord& BamRecord::Name(const std::string& name)
{
    // determine change in memory needed
    // diffNumBytes: pos -> growing, neg -> shrinking
    const size_t numChars = name.size() + 1; // +1 for NULL-term
    const int diffNumBytes = numChars - d_->core.l_qname;
    const int oldLengthData = d_->l_data;
    d_->l_data += diffNumBytes;
    MaybeReallocData();

    // shift trailing data (cigar, seq, qual, tags) as needed
    const uint32_t* oldCigarStart = bam_get_cigar(d_);
    const size_t trailingDataLength = oldLengthData - ((uint8_t*)oldCigarStart - d_->data);
    d_->core.l_qname = numChars;
    uint32_t* newCigarStart = bam_get_cigar(d_);
    memmove(newCigarStart, oldCigarStart, trailingDataLength);

    // fill in new name
    memcpy(d_->data, name.c_str(), numChars);
    return *this;
}

string BamRecord::Qualities(void) const
{
    if (d_->core.l_qseq == 0)
        return string();

    uint8_t* qualData = bam_get_qual(d_);
    if (qualData[0] == 0xff)
        return string();

    string result;
    result.reserve(d_->core.l_qseq);
    for (int i = 0; i < d_->core.l_qseq; ++i)
        result.push_back(qualData[i] + 33);
    return result;
}

/// \cond
std::shared_ptr<bam1_t> BamRecord::RawData(void) const
{
    return d_;
}
/// \endcond

bool BamRecord::RemoveTag(const string& tagName)
{
    if (tagName.size() != 2)
        return false;
    uint8_t* data = bam_aux_get(d_.get(), tagName.c_str());
    if (data == 0)
        return false;
    return bam_aux_del(d_.get(), data) == 0;
}

string BamRecord::Sequence(void) const
{
    string result;
    result.reserve(d_->core.l_qseq);
    static const string DnaLookup = string("=ACMGRSVTWYHKDBN");
    const uint8_t* seqData = bam_get_seq(d_);
    for (int i = 0; i < d_->core.l_qseq; ++i)
        result.append(1, DnaLookup[bam_seqi(seqData, i)]);
    return result;
}

BamRecord& BamRecord::SetSequenceAndQualities(const std::string& sequence,
                                              const std::string& qualities)
{
    // TODO: I'm ok with the assert for now, but how to handle at runtime?
    if (!qualities.empty()) {
        PB_ASSERT_OR_RETURN_VALUE(sequence.size() == qualities.size(), *this);
    }

    return SetSequenceAndQualitiesInternal(sequence.c_str(),
                                           sequence.size(),
                                           qualities.c_str(),
                                           false);
}

BamRecord& BamRecord::SetSequenceAndQualities(const char* sequence,
                                              const size_t sequenceLength,
                                              const char* qualities)
{
    return SetSequenceAndQualitiesInternal(sequence,
                                           sequenceLength,
                                           qualities,
                                           false);
}

BamRecord& BamRecord::SetPreencodedSequenceAndQualities(const char* encodedSequence,
                                                        const size_t rawSequenceLength,
                                                        const char* qualities)
{
    return SetSequenceAndQualitiesInternal(encodedSequence,
                                           rawSequenceLength,
                                           qualities,
                                           true);
}

BamRecord& BamRecord::SetSequenceAndQualitiesInternal(const char* sequence,
                                                      const size_t sequenceLength,
                                                      const char* qualities,
                                                      bool isPreencoded)
{
    // determine change in memory needed
    // diffNumBytes: pos -> growing, neg -> shrinking
    const int encodedSequenceLength = static_cast<int>((sequenceLength+1)/2);
    const int oldSeqAndQualLength = static_cast<int>((d_->core.l_qseq+1)/2) + d_->core.l_qseq; // encoded seq + qual
    const int newSeqAndQualLength = encodedSequenceLength + sequenceLength;                    // encoded seq + qual
    const int diffNumBytes = newSeqAndQualLength - oldSeqAndQualLength;
    const int oldLengthData = d_->l_data;
    d_->l_data += diffNumBytes;
    MaybeReallocData();

    // shift trailing data (tags) as needed
    const uint8_t* oldTagStart = bam_get_aux(d_);
    const size_t trailingDataLength = oldLengthData - ((uint8_t*)oldTagStart - d_->data);
    d_->core.l_qseq = sequenceLength;
    uint8_t* newTagStart = bam_get_aux(d_);
    memmove(newTagStart, oldTagStart, trailingDataLength);

    // fill in new sequence
    uint8_t* pEncodedSequence = bam_get_seq(d_);
    if ( isPreencoded ) {
        memcpy(pEncodedSequence, sequence, encodedSequenceLength);
    } else {
        const char* pRawSequence = sequence;
        uint8_t nucleotideCode;
        bool useHighWord = true;
        for (size_t i = 0; i < sequenceLength; ++i) {
            switch (*pRawSequence) {
                case '=' : nucleotideCode = 0;  break;
                case 'A' : nucleotideCode = 1;  break;
                case 'C' : nucleotideCode = 2;  break;
                case 'M' : nucleotideCode = 3;  break;
                case 'G' : nucleotideCode = 4;  break;
                case 'R' : nucleotideCode = 5;  break;
                case 'S' : nucleotideCode = 6;  break;
                case 'V' : nucleotideCode = 7;  break;
                case 'T' : nucleotideCode = 8;  break;
                case 'W' : nucleotideCode = 9;  break;
                case 'Y' : nucleotideCode = 10; break;
                case 'H' : nucleotideCode = 11; break;
                case 'K' : nucleotideCode = 12; break;
                case 'D' : nucleotideCode = 13; break;
                case 'B' : nucleotideCode = 14; break;
                case 'N' : nucleotideCode = 15; break;
                default :
                    PB_ASSERT_UNREACHABLE; // graceful way to handle?
                    break;
            }

            // pack the nucleotide code
            if (useHighWord) {
                *pEncodedSequence = nucleotideCode << 4;
                useHighWord = false;
            } else {
                *pEncodedSequence |= nucleotideCode;
                ++pEncodedSequence;
                useHighWord = true;
            }
            ++pRawSequence;
        }
    }

    // fill in quality values
    uint8_t* encodedQualities = bam_get_qual(d_);
    if ( (qualities == 0 ) || (::strlen(qualities) == 0) )
        memset(encodedQualities, 0xff, sequenceLength);
    else {
        for (size_t i = 0; i < sequenceLength; ++i)
            encodedQualities[i] = qualities[i] - 33;  // FASTQ ASCII -> int conversion
    }
    return *this;
}

BamRecord& BamRecord::Tags(const TagCollection& tags)
{
    // convert tags to binary
    const vector<uint8_t>& tagData = BamTagCodec::Encode(tags);
    const size_t numBytes = tagData.size();
    const uint8_t* data = tagData.data();

    // determine change in memory needed
    uint8_t* tagStart = bam_get_aux(d_);
    const size_t oldNumBytes = d_->l_data - (tagStart - d_->data);
    const int diffNumBytes = numBytes - oldNumBytes;
    d_->l_data += diffNumBytes;
    MaybeReallocData();
    tagStart = bam_get_aux(d_);

    // fill in new tag data
    memcpy((void*)tagStart, data, numBytes);
    return *this;
}

TagCollection BamRecord::Tags(void) const
{
    const uint8_t* tagDataStart = bam_get_aux(d_);
    const size_t numBytes = d_->l_data - (tagDataStart - d_->data);
    return BamTagCodec::Decode(vector<uint8_t>(tagDataStart, tagDataStart+numBytes));
}

Tag BamRecord::TagValue(const string& tagName) const
{
    if (tagName.size() != 2)
        return Tag();
    uint8_t* data = bam_aux_get(d_.get(), tagName.c_str());
    if (data == 0)
        return Tag();
    return BamTagCodec::FromRawData(data);
}
