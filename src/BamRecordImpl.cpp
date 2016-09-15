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

#include "pbbam/BamRecordImpl.h"
#include "pbbam/BamTagCodec.h"
#include "AssertUtils.h"
#include "BamRecordTags.h"
#include "MemoryUtils.h"
#include <algorithm>
#include <iostream>
#include <utility>
#include <cstdlib>
#include <cstring>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

BamRecordImpl::BamRecordImpl(void)
    : d_(nullptr)
{
    InitializeData();
}

BamRecordImpl::BamRecordImpl(const BamRecordImpl& other)
    : d_(bam_dup1(other.d_.get()), internal::HtslibRecordDeleter())
    , tagOffsets_(other.tagOffsets_)
{ }

BamRecordImpl::BamRecordImpl(BamRecordImpl&& other)
    : d_(nullptr)
    , tagOffsets_(std::move(other.tagOffsets_))
{
    d_.swap(other.d_);
    other.d_.reset();
}

BamRecordImpl& BamRecordImpl::operator=(const BamRecordImpl& other)
{
    if (this != & other) {
        if (d_ == nullptr)
            InitializeData();
        bam_copy1(d_.get(), other.d_.get());
        tagOffsets_ = other.tagOffsets_;
    }
    return *this;
}

BamRecordImpl& BamRecordImpl::operator=(BamRecordImpl&& other)
{
    if (this != & other) {
        d_.swap(other.d_);
        other.d_.reset();

        tagOffsets_ = std::move(other.tagOffsets_);
    }
    return *this;
}

BamRecordImpl::~BamRecordImpl(void) { }

bool BamRecordImpl::AddTag(const string& tagName,
                           const Tag &value)
{
    return AddTag(tagName, value, TagModifier::NONE);
}

bool BamRecordImpl::AddTag(const BamRecordTag tag,
                           const Tag& value)
{
   return AddTag(internal::BamRecordTags::LabelFor(tag),
                 value,
                 TagModifier::NONE);
}

bool BamRecordImpl::AddTag(const string& tagName,
                           const Tag& value,
                           const TagModifier additionalModifier)
{
    if (tagName.size() != 2 || HasTag(tagName))
        return false;
    const bool added = AddTagImpl(tagName, value, additionalModifier);
    if (added)
        UpdateTagMap();
    return added;
}

bool BamRecordImpl::AddTag(const BamRecordTag tag,
                           const Tag& value,
                           const TagModifier additionalModifier)
{
    return AddTag(internal::BamRecordTags::LabelFor(tag),
                  value,
                  additionalModifier);
}

bool BamRecordImpl::AddTagImpl(const string& tagName,
                               const Tag& value,
                               const TagModifier additionalModifier)
{
    const vector<uint8_t> rawData = BamTagCodec::ToRawData(value, additionalModifier);
    if (rawData.empty())
        return false;

    bam_aux_append(d_.get(),
                   tagName.c_str(),
                   BamTagCodec::TagTypeCode(value, additionalModifier),
                   rawData.size(),
                   const_cast<uint8_t*>(rawData.data()));
    return true;
}

Cigar BamRecordImpl::CigarData(void) const
{
    Cigar result;
    result.reserve(d_->core.n_cigar);
    uint32_t* cigarData = bam_get_cigar(d_);
    for (uint32_t i = 0; i < d_->core.n_cigar; ++i) {
        const uint32_t length = bam_cigar_oplen(cigarData[i]);
        const CigarOperationType type = static_cast<CigarOperationType>(bam_cigar_op(cigarData[i]));
        result.push_back(CigarOperation(type, length));
    }

    return result;
}

BamRecordImpl& BamRecordImpl::CigarData(const Cigar& cigar)
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
        cigarDataStart[i] = bam_cigar_gen(cigarOp.Length(), static_cast<int>(cigarOp.Type()));
    }

    return *this;
}

BamRecordImpl& BamRecordImpl::CigarData(const std::string& cigarString)
{
    return CigarData(Cigar::FromStdString(cigarString));
}

bool BamRecordImpl::EditTag(const string& tagName,
                            const Tag& newValue)
{
    return EditTag(tagName, newValue, TagModifier::NONE);
}

bool BamRecordImpl::EditTag(const BamRecordTag tag,
                            const Tag& newValue)
{
    return EditTag(internal::BamRecordTags::LabelFor(tag),
                   newValue,
                   TagModifier::NONE);
}

bool BamRecordImpl::EditTag(const string& tagName,
                            const Tag& newValue,
                            const TagModifier additionalModifier)
{
    // try remove old value (with delayed tag map update)
    const bool removed = RemoveTagImpl(tagName);
    if (!removed)
        return false;

    // if old value removed, add new value
    const bool added = AddTagImpl(tagName, newValue, additionalModifier);
    if (added)
        UpdateTagMap();
    return added;
}

bool BamRecordImpl::EditTag(const BamRecordTag tag,
                            const Tag& newValue,
                            const TagModifier additionalModifier)
{
    return EditTag(internal::BamRecordTags::LabelFor(tag),
                   newValue,
                   additionalModifier);
}

BamRecordImpl BamRecordImpl::FromRawData(const PBBAM_SHARED_PTR<bam1_t>& rawData)
{
    BamRecordImpl result;
    bam_copy1(result.d_.get(), rawData.get());
    return result;
}

bool BamRecordImpl::HasTag(const string& tagName) const
{
    if (tagName.size() != 2)
        return false;
    return TagOffset(tagName) != -1;

    // 27635
//    return bam_aux_get(d_.get(), tagName.c_str()) != 0;
}

bool BamRecordImpl::HasTag(const BamRecordTag tag) const
{
    return HasTag(internal::BamRecordTags::LabelFor(tag));
}

void BamRecordImpl::InitializeData(void)
{
    d_.reset(bam_init1(), internal::HtslibRecordDeleter());
    d_->data = (uint8_t*)(calloc(0x800, sizeof(uint8_t)));   // maybe make this value tune-able later?

    // init unmapped
    Position(PacBio::BAM::UnmappedPosition);
    MatePosition(PacBio::BAM::UnmappedPosition);
    ReferenceId(-1);
    MateReferenceId(-1);
    SetMapped(false);
    MapQuality(255);

    // initialized with NULL term for qname
    d_->core.l_qname = 1;
    d_->l_data = 1;
    d_->m_data = 0x800;
}

void BamRecordImpl::MaybeReallocData(void)
{
    // about to grow data contents to l_data size, but m_data is our current max.
    // so we may need to grow. if so, use kroundup to double to next power of 2
    if (d_->m_data < d_->l_data) {
        d_->m_data = d_->l_data;
        kroundup32(d_->m_data);
        d_->data = static_cast<uint8_t*>(realloc(d_->data, d_->m_data));
    }
}

string BamRecordImpl::Name(void) const
{
    return string(bam_get_qname(d_));
}

BamRecordImpl& BamRecordImpl::Name(const std::string& name)
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

QualityValues BamRecordImpl::Qualities(void) const
{
    if (d_->core.l_qseq == 0)
        return QualityValues();

    uint8_t* qualData = bam_get_qual(d_);
    if (qualData[0] == 0xff)
        return QualityValues();

    const size_t numQuals = d_->core.l_qseq;
    QualityValues result;
    result.reserve(numQuals);
    for (size_t i = 0; i < numQuals; ++i)
        result.push_back(QualityValue(qualData[i]));
    return result;
}

bool BamRecordImpl::RemoveTag(const string& tagName)
{
    const bool removed = RemoveTagImpl(tagName);
    if (removed)
        UpdateTagMap();
    return removed;
}

bool BamRecordImpl::RemoveTag(const BamRecordTag tag)
{
    return RemoveTag(internal::BamRecordTags::LabelFor(tag));
}

bool BamRecordImpl::RemoveTagImpl(const string &tagName)
{
    if (tagName.size() != 2)
        return false;
    uint8_t* data = bam_aux_get(d_.get(), tagName.c_str());
    if (data == 0)
        return false;
    const bool ok = bam_aux_del(d_.get(), data) == 0;
    return ok;
}

string BamRecordImpl::Sequence(void) const
{
    string result;
    result.reserve(d_->core.l_qseq);
    static const string DnaLookup = string("=ACMGRSVTWYHKDBN");
    const uint8_t* seqData = bam_get_seq(d_);
    for (int i = 0; i < d_->core.l_qseq; ++i)
        result.append(1, DnaLookup[bam_seqi(seqData, i)]);
    return result;
}

size_t BamRecordImpl::SequenceLength(void) const
{ return d_->core.l_qseq; }

BamRecordImpl& BamRecordImpl::SetSequenceAndQualities(const std::string& sequence,
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

BamRecordImpl& BamRecordImpl::SetSequenceAndQualities(const char* sequence,
                                                      const size_t sequenceLength,
                                                      const char* qualities)
{
    return SetSequenceAndQualitiesInternal(sequence,
                                           sequenceLength,
                                           qualities,
                                           false);
}

BamRecordImpl& BamRecordImpl::SetPreencodedSequenceAndQualities(const char* encodedSequence,
                                                                const size_t rawSequenceLength,
                                                                const char* qualities)
{
    return SetSequenceAndQualitiesInternal(encodedSequence,
                                           rawSequenceLength,
                                           qualities,
                                           true);
}

BamRecordImpl& BamRecordImpl::SetSequenceAndQualitiesInternal(const char* sequence,
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
    if (isPreencoded) {
        memcpy(pEncodedSequence, sequence, encodedSequenceLength);
    } else {
        memset(pEncodedSequence, 0, encodedSequenceLength);
        for (size_t i = 0; i < sequenceLength; ++i)
            pEncodedSequence[i>>1] |= seq_nt16_table[(int)sequence[i]] << ((~i&1)<<2);
    }

    // fill in quality values
    uint8_t* encodedQualities = bam_get_qual(d_);
    if ( (qualities == 0 ) || (strlen(qualities) == 0) )
        memset(encodedQualities, 0xff, sequenceLength);
    else {
        for (size_t i = 0; i < sequenceLength; ++i)
            encodedQualities[i] = qualities[i] - 33;  // FASTQ ASCII -> int conversion
    }
    return *this;
}

int BamRecordImpl::TagOffset(const string& tagName) const
{
    if (tagName.size() != 2)
        throw std::runtime_error("invalid tag name size");

    if (tagOffsets_.empty())
        UpdateTagMap();

    const uint16_t tagCode = (static_cast<uint8_t>(tagName.at(0)) << 8) | static_cast<uint8_t>(tagName.at(1));
    const auto found = tagOffsets_.find(tagCode);
    return (found != tagOffsets_.cend() ? found->second : -1);
}

BamRecordImpl& BamRecordImpl::Tags(const TagCollection& tags)
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

    // update tag info
    UpdateTagMap();
    return *this;
}

TagCollection BamRecordImpl::Tags(void) const
{
    const uint8_t* tagDataStart = bam_get_aux(d_);
    const size_t numBytes = d_->l_data - (tagDataStart - d_->data);
    return BamTagCodec::Decode(vector<uint8_t>(tagDataStart, tagDataStart+numBytes));
}

Tag BamRecordImpl::TagValue(const string& tagName) const
{
    if (tagName.size() != 2)
        return Tag();

    const int offset = TagOffset(tagName);
    if (offset == -1)
        return Tag();

    bam1_t* b = d_.get();
    assert(bam_get_aux(b));
    uint8_t* tagData = bam_get_aux(b) + offset;
    if (offset >= b->l_data)
        return Tag();

    // skip tag name
    return BamTagCodec::FromRawData(tagData);
}

Tag BamRecordImpl::TagValue(const BamRecordTag tag) const
{
    return TagValue(internal::BamRecordTags::LabelFor(tag));
}

void BamRecordImpl::UpdateTagMap(void) const
{
    // clear out offsets, leave map structure basically intact
    auto tagIter = tagOffsets_.begin();
    auto tagEnd  = tagOffsets_.end();
    for ( ; tagIter != tagEnd; ++tagIter )
        tagIter->second = -1;

    const uint8_t* tagStart = bam_get_aux(d_);
    if (tagStart == 0)
        return;
    const ptrdiff_t numBytes = d_->l_data - (tagStart - d_->data);

    // NOTE: using a 16-bit 'code' for tag name here instead of string, to avoid
    // a lot of string constructions & comparisons. All valid tags will be 2 chars
    // anyway, so this should be a nice lookup mechanism.
    //
    uint16_t tagNameCode;
    int64_t i = 0;
    while(i < numBytes) {

        // store (tag name code -> start offset into tag data)
        tagNameCode = static_cast<char>(tagStart[i]) << 8 | static_cast<char>(tagStart[i+1]);
        i += 2;
        tagOffsets_[tagNameCode] = i;

        // skip tag contents
        const char tagType = static_cast<char>(tagStart[i++]);
        switch (tagType) {
            case 'A' :
            case 'a' :
            case 'c' :
            case 'C' :
            {
                i += 1;
                break;
            }
            case 's' :
            case 'S' :
            {
                i += 2;
                break;
            }
            case 'i' :
            case 'I' :
            case 'f' :
            {
                i += 4;
                break;
            }

            case 'Z' :
            case 'H' :
            {
                // null-terminated string
                i += strlen((const char*)&tagStart[i]) + 1;
                break;
            }

            case 'B' :
            {
                const char subTagType = tagStart[i++];
                size_t elementSize = 0;
                switch (subTagType) {
                    case 'c' :
                    case 'C' : elementSize = 1; break;
                    case 's' :
                    case 'S' : elementSize = 2; break;
                    case 'i' :
                    case 'I' :
                    case 'f' : elementSize = 4; break;

                    // unknown subTagType
                    default:
                        PB_ASSERT_OR_RETURN(false);
                }

                uint32_t numElements = 0;
                memcpy(&numElements, &tagStart[i], sizeof(uint32_t));
                i += (4 + (elementSize * numElements));
                break;
            }

            // unknown tagType
            default:
                PB_ASSERT_OR_RETURN(false);
        }
    }
}
