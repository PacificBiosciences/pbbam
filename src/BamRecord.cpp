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
/// \file BamRecord.cpp
/// \brief Implements the BamRecord & BamRecordView classes.
//
// Author: Derek Barnett

#include "pbbam/BamRecord.h"
#include "pbbam/virtual/VirtualRegionTypeMap.h"
#include "pbbam/ZmwTypeMap.h"
#include "AssertUtils.h"
#include "MemoryUtils.h"
#include "SequenceUtils.h"
#include <boost/numeric/conversion/cast.hpp>
#include <htslib/sam.h>
#include <iostream>
#include <stdexcept>

using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

// BAM record tag names
static const string tagName_alternative_labelQV     = "pv";
static const string tagName_alternative_labelTag    = "pt";
static const string tagName_barcodes                = "bc";
static const string tagName_barcode_quality         = "bq";
static const string tagName_contextFlags            = "cx";
static const string tagName_holeNumber              = "zm";
static const string tagName_deletionQV              = "dq";
static const string tagName_deletionTag             = "dt";
static const string tagName_insertionQV             = "iq";
static const string tagName_ipd                     = "ip";
static const string tagName_labelQV                 = "pq";
static const string tagName_mergeQV                 = "mq";
static const string tagName_numPasses               = "np";
static const string tagName_pkmean                  = "pa";
static const string tagName_pkmid                   = "pm";
static const string tagName_pkmean2                 = "ps";
static const string tagName_pkmid2                  = "pi";
static const string tagName_pre_pulse_frames        = "pd";
static const string tagName_pulse_call              = "pc";
static const string tagName_pulse_call_width        = "px";
static const string tagName_pulseMergeQV            = "pg";
static const string tagName_pulseWidth              = "pw";
static const string tagName_queryStart              = "qs";
static const string tagName_queryEnd                = "qe";
static const string tagName_readAccuracy            = "rq";
static const string tagName_readGroup               = "RG";
static const string tagName_scrap_region_type       = "sc";
static const string tagName_scrap_zmw_type          = "sz";
static const string tagName_snr                     = "sn";
static const string tagName_startFrame              = "sf";
static const string tagName_substitutionQV          = "sq";
static const string tagName_substitutionTag         = "st";

// faux (helper) tag names
static const string tagName_QUAL = "QUAL";
static const string tagName_SEQ  = "SEQ";

// record type names
static const string recordTypeName_ZMW        = "ZMW";
static const string recordTypeName_Polymerase = "POLYMERASE";
static const string recordTypeName_HqRegion   = "HQREGION";
static const string recordTypeName_Subread    = "SUBREAD";
static const string recordTypeName_CCS        = "CCS";
static const string recordTypeName_Scrap      = "SCRAP";
static const string recordTypeName_Unknown    = "UNKNOWN";

static
int32_t HoleNumberFromName(const string& fullName)
{
    const auto mainTokens = Split(fullName, '/');
    if (mainTokens.size() != 3)
        throw std::runtime_error("malformed record name");
    return stoi(mainTokens.at(1));
}

static
Position QueryEndFromName(const string& fullName)
{
    const auto mainTokens = Split(fullName, '/');
    if (mainTokens.size() != 3)
        throw std::runtime_error("malformed record name");
    const auto queryTokens = Split(mainTokens.at(2), '_');
    if (queryTokens.size() != 2)
        throw std::runtime_error("malformed record name");
    return stoi(queryTokens.at(1));
}

static
Position QueryStartFromName(const string& fullName)
{
    const auto mainTokens = Split(fullName, '/');
    if (mainTokens.size() != 3)
        throw std::runtime_error("malformed record name");
    const auto queryTokens = Split(mainTokens.at(2), '_');
    if (queryTokens.size() != 2)
        throw std::runtime_error("malformed record name");
    return stoi(queryTokens.at(0));
}

static
BamRecordImpl* CreateOrEdit(const string& tagName,
                            const Tag& value,
                            BamRecordImpl* impl)
{
    if (impl->HasTag(tagName))
        impl->EditTag(tagName, value);
    else
        impl->AddTag(tagName, value);
    return impl;
}

static
pair<int32_t, int32_t> AlignedOffsets(const BamRecord& record,
                                      const int seqLength)
{
    int32_t startOffset = 0;
    int32_t endOffset = seqLength;

    PBBAM_SHARED_PTR<bam1_t> b = internal::BamRecordMemory::GetRawData(record);
    uint32_t* cigarData = bam_get_cigar(b.get());
    const size_t numCigarOps = b->core.n_cigar;
    if (numCigarOps > 0) {

        // start offset
        for (size_t i = 0; i < numCigarOps; ++i) {
            const CigarOperationType type = static_cast<CigarOperationType>(bam_cigar_op(cigarData[i]));
            if (type == CigarOperationType::HARD_CLIP) {
                if (startOffset != 0 && startOffset != seqLength) {
                    startOffset = -1;
                    break;
                }
            }
            else if (type == CigarOperationType::SOFT_CLIP)
                startOffset += bam_cigar_oplen(cigarData[i]);
            else
                break;
        }

        // end offset
        for (int i = numCigarOps-1; i >= 0; --i) {
            const CigarOperationType type = static_cast<CigarOperationType>(bam_cigar_op(cigarData[i]));
            if (type == CigarOperationType::HARD_CLIP) {
                if (endOffset != 0 && endOffset != seqLength) {
                    endOffset = -1;
                    break;
                }
            }
            else if (type == CigarOperationType::SOFT_CLIP)
                endOffset -= bam_cigar_oplen(cigarData[i]);
            else
                break;

        }

        if (endOffset == 0)
            endOffset = seqLength;
    }
    return std::make_pair(startOffset, endOffset);
}

template<typename T>
T Clip(const T& input, const size_t pos, const size_t len)
{
    if (input.empty())
        return T();
    return T{ input.cbegin() + pos,
              input.cbegin() + pos + len };
}

static
void ClipAndGapifyBases(const BamRecordImpl& impl,
                        const bool aligned,
                        const bool exciseSoftClips,
                        string* seq)
{
    assert(seq);
    if (impl.IsMapped() && (aligned || exciseSoftClips)) {

        size_t seqIndex = 0;
        const auto cigar = impl.CigarData();
        auto cigarIter = cigar.cbegin();
        auto cigarEnd  = cigar.cend();
        for (; cigarIter != cigarEnd; ++cigarIter) {
            const auto op = (*cigarIter);
            const auto type = op.Type();

            // do nothing for hard clips
            if (type != CigarOperationType::HARD_CLIP) {
                const auto opLength = op.Length();

                // maybe remove soft clips
                if (type == CigarOperationType::SOFT_CLIP && exciseSoftClips)
                    seq->erase(seqIndex, opLength);

                // for non-clipping operations
                else {

                    // maybe add gaps/padding
                    if (aligned) {
                        if (type == CigarOperationType::DELETION) {
                            seq->reserve(seq->size() + opLength);
                            seq->insert(seqIndex, opLength, '-');
                        }
                        else if (type == CigarOperationType::PADDING) {
                            seq->reserve(seq->size() + opLength);
                            seq->insert(seqIndex, opLength, '*');
                        }
                    }

                    // update index
                    seqIndex += opLength;
                }
            }
        }
    }
}

static
void ClipAndGapifyFrames(const BamRecordImpl& impl,
                         const bool aligned,
                         const bool exciseSoftClips,
                         Frames* frames)
{
    assert(frames);

    if (impl.IsMapped() && (aligned || exciseSoftClips)) {

        auto data = std::move(frames->Data()); // we're going to put it back
        size_t frameIndex = 0;
        const auto cigar = impl.CigarData();
        auto cigarIter = cigar.cbegin();
        auto cigarEnd  = cigar.cend();
        for (; cigarIter != cigarEnd; ++cigarIter) {
            const auto op = (*cigarIter);
            const auto type = op.Type();

            // do nothing for hard clips
            if (type != CigarOperationType::HARD_CLIP) {
                const auto opLength = op.Length();

                // maybe remove soft clips
                if (type == CigarOperationType::SOFT_CLIP && exciseSoftClips) {
                    data.erase(data.begin() + frameIndex,
                               data.begin() + frameIndex + opLength);
                }

                // for non-clipping operations
                else {

                    // maybe add gaps/padding
                    if (aligned) {
                        if (type == CigarOperationType::DELETION ||
                            type == CigarOperationType::PADDING)
                        {
                            data.reserve(data.size() + opLength);
                            data.insert(data.begin() + frameIndex, opLength, 0);
                        }
                    }

                    // update index
                    frameIndex += opLength;
                }
            }
        }
        frames->Data(data);
    }
}


static
void ClipAndGapifyQualities(const BamRecordImpl& impl,
                            const bool aligned,
                            const bool exciseSoftClips,
                            QualityValues* quals)
{
    assert(quals);

    if (impl.IsMapped() && (aligned || exciseSoftClips)) {

        size_t qualIndex = 0;
        const auto cigar = impl.CigarData();
        auto cigarIter = cigar.cbegin();
        auto cigarEnd  = cigar.cend();
        for (; cigarIter != cigarEnd; ++cigarIter) {
            const auto op = (*cigarIter);
            const auto type = op.Type();

            // do nothing for hard clips
            if (type != CigarOperationType::HARD_CLIP) {
                const auto opLength = op.Length();

                // maybe remove soft clips
                if (type == CigarOperationType::SOFT_CLIP && exciseSoftClips)
                    quals->erase(quals->begin() + qualIndex, quals->begin() + qualIndex + opLength);

                // for non-clipping operations
                else {

                    // maybe add gaps/padding
                    if (aligned) {
                        if (type == CigarOperationType::DELETION || type == CigarOperationType::PADDING) {
                            quals->reserve(quals->size() + opLength);
                            quals->insert(quals->begin() + qualIndex, opLength, QualityValue(0));
                        }
                    }

                    // update index
                    qualIndex += opLength;
                }
            }
        }
    }
}

static
RecordType NameToType(const string& name)
{
    if (name == recordTypeName_Subread)
        return RecordType::SUBREAD;
    if (name == recordTypeName_ZMW || name == recordTypeName_Polymerase)
        return RecordType::ZMW;
    if (name == recordTypeName_HqRegion)
        return RecordType::HQREGION;
    if (name == recordTypeName_CCS)
        return RecordType::CCS;
    if (name == recordTypeName_Scrap)
        return RecordType::SCRAP;
    return RecordType::UNKNOWN;
}

static
void OrientBasesAsRequested(string* bases,
                            Orientation current,
                            Orientation requested,
                            bool isReverseStrand,
                            bool isPulse)
{
    assert(bases);
    if (current != requested && isReverseStrand) {
        if (isPulse)
            internal::ReverseComplementCaseSens(*bases);
        else
            internal::ReverseComplement(*bases);
    }
}

template<typename Container> inline
void OrientTagDataAsRequested(Container* data,
                              Orientation current,
                              Orientation requested,
                              bool isReverseStrand)
{
    assert(data);
    if (current != requested && isReverseStrand)
        std::reverse(data->begin(), data->end());
}

static inline
bool ConsumesQuery(const CigarOperationType type)
{ return (bam_cigar_type(static_cast<int>(type)) & 0x1) != 0; }

static inline
bool ConsumesReference(const CigarOperationType type)
{ return (bam_cigar_type(static_cast<int>(type)) & 0x2) != 0; }

} // namespace internal
} // namespace BAM
} // namespace PacBio

const float BamRecord::photonFactor = 10.0;

BamRecord::BamRecord(void)
    : alignedStart_(PacBio::BAM::UnmappedPosition)
    , alignedEnd_(PacBio::BAM::UnmappedPosition)
{ }

BamRecord::BamRecord(const BamHeader& header)
    : header_(header)
    , alignedStart_(PacBio::BAM::UnmappedPosition)
    , alignedEnd_(PacBio::BAM::UnmappedPosition)
{ }

BamRecord::BamRecord(const BamRecordImpl& impl)
    : impl_(impl)
    , alignedStart_(PacBio::BAM::UnmappedPosition)
    , alignedEnd_(PacBio::BAM::UnmappedPosition)
{ }

BamRecord::BamRecord(BamRecordImpl&& impl)
    : impl_(std::move(impl))
    , alignedStart_(PacBio::BAM::UnmappedPosition)
    , alignedEnd_(PacBio::BAM::UnmappedPosition)
{ }

BamRecord::BamRecord(const BamRecord& other)
    : impl_(other.impl_)
    , header_(other.header_)
    , alignedStart_(other.alignedStart_)
    , alignedEnd_(other.alignedEnd_)
{ }

BamRecord::BamRecord(BamRecord&& other)
    : impl_(std::move(other.impl_))
    , header_(std::move(other.header_))
    , alignedStart_(std::move(other.alignedStart_))
    , alignedEnd_(std::move(other.alignedEnd_))
{ }

BamRecord& BamRecord::operator=(const BamRecord& other)
{
    impl_ = other.impl_;
    header_ = other.header_;
    alignedStart_ = other.alignedStart_;
    alignedEnd_ = other.alignedEnd_;
    return *this;
}

BamRecord& BamRecord::operator=(BamRecord&& other)
{
    impl_ = std::move(other.impl_);
    header_ = std::move(other.header_);
    alignedStart_ = std::move(other.alignedStart_);
    alignedEnd_ = std::move(other.alignedEnd_);
    return *this;
}

BamRecord::~BamRecord(void) { }

Position BamRecord::AlignedEnd(void) const
{
    if (alignedEnd_ == PacBio::BAM::UnmappedPosition)
        CalculateAlignedPositions();
    return alignedEnd_;
}

Position BamRecord::AlignedStart(void) const
{
    if (alignedStart_ == PacBio::BAM::UnmappedPosition)
        CalculateAlignedPositions();
    return alignedStart_;
}

Strand BamRecord::AlignedStrand(void) const
{ return impl_.IsReverseStrand() ? Strand::REVERSE : Strand::FORWARD; }

QualityValues BamRecord::AltLabelQV(Orientation orientation) const
{
    return FetchQualities(internal::tagName_alternative_labelQV,
                          orientation);
}

BamRecord& BamRecord::AltLabelQV(const QualityValues& altLabelQVs)
{
    internal::CreateOrEdit(internal::tagName_alternative_labelQV,
                           altLabelQVs.Fastq(), &impl_);
    return *this;
}

std::string BamRecord::AltLabelTag(Orientation orientation) const
{
    return FetchBases(internal::tagName_alternative_labelTag,
                      orientation);
}

BamRecord& BamRecord::AltLabelTag(const std::string& tags)
{
    internal::CreateOrEdit(internal::tagName_alternative_labelTag, tags, &impl_);
    return *this;
}

int16_t BamRecord::BarcodeForward(void) const
{ return Barcodes().first; }

int16_t BamRecord::BarcodeReverse(void) const
{ return Barcodes().second; }

uint8_t BamRecord::BarcodeQuality(void) const
{
    const auto bq = impl_.TagValue(internal::tagName_barcode_quality);
    if (bq.IsNull())
        return 0; // ?? "missing" value for tags ?? should we consider boost::optional<T> for these kind of guys ??
    return bq.ToUInt8();
}

BamRecord& BamRecord::BarcodeQuality(const uint8_t quality)
{
    internal::CreateOrEdit(internal::tagName_barcode_quality, quality, &impl_);
    return *this;
}

std::pair<int16_t,int16_t> BamRecord::Barcodes(void) const
{
    const Tag& bc = impl_.TagValue(internal::tagName_barcodes);
    if (bc.IsNull())
        throw std::runtime_error("barcode tag (bc) was requested but is missing");

    // NOTE: barcodes are still stored, per the spec, as uint16, even though
    // we're now using them as int16_t in the API (bug 31511)
    //
    if (!bc.IsUInt16Array())
        throw std::runtime_error("barcode tag (bc) is malformed: should be a uint16_t array of size==2.");
    const auto bcArray = bc.ToUInt16Array();
    if (bcArray.size() != 2)
        throw std::runtime_error("barcode tag (bc) is malformed: should be a uint16_t array of size==2.");

    return std::make_pair(boost::numeric_cast<int16_t>(bcArray[0]),
                          boost::numeric_cast<int16_t>(bcArray[1]));
}

BamRecord& BamRecord::Barcodes(const std::pair<int16_t,int16_t>& barcodeIds)
{
    const std::vector<uint16_t> data =
    {
        boost::numeric_cast<uint16_t>(barcodeIds.first),
        boost::numeric_cast<uint16_t>(barcodeIds.second)
    };
    internal::CreateOrEdit(internal::tagName_barcodes, data, &impl_);
    return *this;
}

void BamRecord::CalculateAlignedPositions(void) const
{
    // reset
    ResetCachedPositions();

    // skip if unmapped, or has no queryStart/End
    if (!impl_.IsMapped())
        return;

    // get the query start/end
    const size_t seqLength = impl_.SequenceLength();
    const RecordType type  = Type();
    const Position qStart  = (type == RecordType::CCS) ? Position(0) : QueryStart();
    const Position qEnd    = (type == RecordType::CCS) ? Position(seqLength) : QueryEnd();
    
    if (qStart == PacBio::BAM::UnmappedPosition || qEnd == PacBio::BAM::UnmappedPosition)
        return;

    // determine clipped end ranges
    const std::pair<int32_t, int32_t> alignedOffsets = internal::AlignedOffsets(*this, seqLength);
    const int32_t startOffset = alignedOffsets.first;
    const int32_t endOffset = alignedOffsets.second;
    if (endOffset == -1 || startOffset == -1)
        return; // TODO: handle error more??

    // store aligned positions (polymerase read coordinates)
    if (impl_.IsReverseStrand()) {
        alignedStart_ = qStart + (seqLength - endOffset);
        alignedEnd_   = qEnd - startOffset;
    }
    else {
        alignedStart_ = qStart + startOffset;
        alignedEnd_   = qEnd - (seqLength - endOffset);
    }
}

Cigar BamRecord::CigarData(bool exciseAllClips) const
{
    auto isClippingOp = [](const CigarOperation& op)
    {
        const auto type = op.Type();
        return type == CigarOperationType::SOFT_CLIP ||
               type == CigarOperationType::HARD_CLIP;
    };

    auto cigar = impl_.CigarData();
    if (exciseAllClips) {
        cigar.erase(std::remove_if(cigar.begin(),
                                   cigar.end(),
                                   isClippingOp),
                    cigar.end());
    }
    return cigar;
}

BamRecord& BamRecord::Clip(const ClipType clipType,
                           const Position start,
                           const Position end)
{
    switch (clipType)
    {
        case ClipType::CLIP_NONE         : return *this;
        case ClipType::CLIP_TO_QUERY     : return ClipToQuery(start, end);
        case ClipType::CLIP_TO_REFERENCE : return ClipToReference(start, end);
        default:
            throw std::runtime_error("unsupported clip type requested");
    }
}

void BamRecord::ClipFields(const size_t clipFrom,
                           const size_t clipLength)
{
    const bool isForwardStrand = (AlignedStrand() == Strand::FORWARD);

    // clip seq, quals
    string sequence = std::move(internal::Clip(Sequence(Orientation::NATIVE), clipFrom, clipLength));
    QualityValues qualities = std::move(internal::Clip(Qualities(Orientation::NATIVE), clipFrom, clipLength));
    if (!isForwardStrand) {
        internal::ReverseComplement(sequence);
        internal::Reverse(qualities);
    }
    impl_.SetSequenceAndQualities(sequence, qualities.Fastq());

    // update BAM tags
    TagCollection tags = impl_.Tags();
    if (HasDeletionQV())
        tags[internal::tagName_deletionQV]          = internal::Clip(DeletionQV(Orientation::NATIVE), clipFrom, clipLength).Fastq();
    if (HasInsertionQV())
        tags[internal::tagName_insertionQV]         = internal::Clip(InsertionQV(Orientation::NATIVE), clipFrom, clipLength).Fastq();
    if (HasMergeQV())
        tags[internal::tagName_mergeQV]             = internal::Clip(MergeQV(Orientation::NATIVE), clipFrom, clipLength).Fastq();
    if (HasSubstitutionQV())
        tags[internal::tagName_substitutionQV]      = internal::Clip(SubstitutionQV(Orientation::NATIVE), clipFrom, clipLength).Fastq();
    if (HasIPD())
        tags[internal::tagName_ipd]                 = internal::Clip(IPD(Orientation::NATIVE).Data(), clipFrom, clipLength);
    if (HasPulseWidth())
        tags[internal::tagName_pulseWidth]          = internal::Clip(PulseWidth(Orientation::NATIVE).Data(), clipFrom, clipLength);
    if (HasDeletionTag())
        tags[internal::tagName_deletionTag]         = internal::Clip(DeletionTag(Orientation::NATIVE), clipFrom, clipLength);
    if (HasSubstitutionTag())
        tags[internal::tagName_substitutionTag]     = internal::Clip(SubstitutionTag(Orientation::NATIVE), clipFrom, clipLength);

    // need to implement clipping on pulse tags etc (bug 31633)
    if (HasAltLabelQV())
        tags[internal::tagName_alternative_labelQV] = AltLabelQV(Orientation::NATIVE).Fastq();
    if (HasLabelQV())
        tags[internal::tagName_labelQV]             = LabelQV(Orientation::NATIVE).Fastq();
    if (HasPulseMergeQV())
        tags[internal::tagName_pulseMergeQV]        = PulseMergeQV(Orientation::NATIVE).Fastq();
    if (HasAltLabelTag())
        tags[internal::tagName_alternative_labelTag]= AltLabelTag(Orientation::NATIVE);
    if (HasPulseCall())
        tags[internal::tagName_pulse_call]          = PulseCall(Orientation::NATIVE);
    if (HasPkmean())
        tags[internal::tagName_pkmean]              = EncodePhotons(Pkmean(Orientation::NATIVE));
    if (HasPkmid())
        tags[internal::tagName_pkmid]               = EncodePhotons(Pkmid(Orientation::NATIVE));
    if (HasPkmean2())
        tags[internal::tagName_pkmean2]             = EncodePhotons(Pkmean2(Orientation::NATIVE));
    if (HasPkmid2())
        tags[internal::tagName_pkmid2]              = EncodePhotons(Pkmid2(Orientation::NATIVE));
    if (HasPrePulseFrames())
        tags[internal::tagName_pre_pulse_frames]    = PrePulseFrames(Orientation::NATIVE).Data();
    if (HasPulseCallWidth())
        tags[internal::tagName_pulse_call_width]    = PulseCallWidth(Orientation::NATIVE).Data();
    if (HasStartFrame())
        tags[internal::tagName_startFrame]          = StartFrame(Orientation::NATIVE);

    impl_.Tags(tags);
}

BamRecord& BamRecord::ClipToQuery(const Position start,
                                  const Position end)
{
    // cache original coords, skip out if clip not needed
    const Position origQStart = QueryStart();
    const Position origQEnd   = QueryEnd();
    if (start <= origQStart && end >= origQEnd)
        return *this;

    // determine new offsets into data
    const size_t startOffset = start - origQStart;
    const size_t endOffset   = origQEnd - end;

    // maybe update CIGAR & aligned position
    if (IsMapped()) {

        // fetch a 'working copy' of CIGAR data
        Cigar cigar = std::move(impl_.CigarData());

        // clip leading CIGAR ops
        size_t referencePositionOffset = 0;
        size_t remaining = startOffset;
        while (remaining > 0 && !cigar.empty()) {
            CigarOperation& firstOp = cigar.front();
            const size_t firstOpLength = firstOp.Length();
            const bool consumesQuery = internal::ConsumesQuery(firstOp.Type());
            const bool consumesRef = internal::ConsumesReference(firstOp.Type());

            // if (!consumesQuery)
            //    just pop (e.g. deletion) ?
            // else {
            //    check bounds, like clip to reference ?
            // }

            // CIGAR op ends at or before clip
            if (firstOpLength <= remaining) {
                cigar.erase(cigar.begin());
                if (consumesQuery)
                    remaining -= firstOpLength;
                if (consumesRef)
                    referencePositionOffset += firstOpLength;
            }

            // CIGAR op straddles clip
            else {
                firstOp.Length(firstOpLength - remaining);
                if (consumesRef)
                    referencePositionOffset += remaining;
                remaining = 0;
            }
        }

        // clip trailing CIGAR ops
        remaining = endOffset;
        while (remaining > 0 && !cigar.empty()) {
            CigarOperation& lastOp = cigar.back();
            const size_t lastOpLength = lastOp.Length();
            const bool consumesQuery = internal::ConsumesQuery(lastOp.Type());

            // CIGAR op ends at or after clip
            if (lastOpLength <= remaining) {
                cigar.pop_back();
                if (consumesQuery)
                    remaining -= lastOpLength;
            }

            // CIGAR op straddles clip
            else {
                lastOp.Length(lastOpLength - remaining);
                remaining = 0;
            }
        }

        // update CIGAR & position
        impl_.CigarData(cigar);
        const Position origPosition = impl_.Position();
        impl_.Position(origPosition + referencePositionOffset);
    }

    // clip SEQ, QUAL, & tags
    const size_t clipFrom   = startOffset;
    const size_t clipLength = (end - start);
    ClipFields(clipFrom, clipLength);

    // update query start/end
    // TODO: update name to reflect new QS/QE ???
    internal::CreateOrEdit(internal::tagName_queryStart, start, &impl_);
    internal::CreateOrEdit(internal::tagName_queryEnd,   end,   &impl_);
//    UpdateName();

    // reset any cached aligned start/end
    ResetCachedPositions();
    return *this;
}

BamRecord& BamRecord::ClipToReference(const Position start,
                                      const Position end)
{
    // skip if not mapped, clipping to reference doesn't make sense
    // or should we even consider throwing here?
    if (!IsMapped())
        return *this;

    const bool isForwardStrand = (AlignedStrand() == Strand::FORWARD);
    return (isForwardStrand ? ClipToReferenceForward(start, end)
                            : ClipToReferenceReverse(start, end));
}

BamRecord& BamRecord::ClipToReferenceForward(const PacBio::BAM::Position start,
                                             const PacBio::BAM::Position end)
{
    assert(IsMapped());
    assert(AlignedStrand() == Strand::FORWARD);

    // cache original coords
    const Position origQStart = QueryStart();
    const Position origQEnd   = QueryEnd();
    const Position origTStart = ReferenceStart();
    const Position origTEnd   = ReferenceEnd();
    assert(AlignedStart() >= origQStart);
    assert(AlignedEnd()   <= origQEnd);

    // skip if already within requested clip range
    if (start <= origTStart && end >= origTEnd)
        return *this;

    const Position newTStart = std::max(origTStart, start);
    const Position newTEnd   = std::min(origTEnd, end);

    // fetch a 'working copy' of CIGAR data
    Cigar cigar = std::move(impl_.CigarData());

    // we're going to skip query sequence outside aligned region
    size_t queryPosRemovedFront = 0;
    size_t queryPosRemovedBack  = 0;

    // ------------------------
    // clip leading CIGAR ops
    // ------------------------

    size_t remaining = newTStart - origTStart;
    while (remaining > 0 && !cigar.empty()) {
        CigarOperation& firstOp = cigar.front();
        const size_t firstOpLength = firstOp.Length();
        const bool consumesQuery = internal::ConsumesQuery(firstOp.Type());
        const bool consumesRef = internal::ConsumesReference(firstOp.Type());

        if (!consumesRef) {

            // e.g. softclip - just pop it completely
            cigar.erase(cigar.begin());
            if (consumesQuery)
                queryPosRemovedFront += firstOpLength;

        } else {
            assert(consumesRef);

            // CIGAR ends at or before clip
            if (firstOpLength <= remaining) {
                cigar.erase(cigar.begin());
                if (consumesQuery)
                    queryPosRemovedFront += firstOpLength;
                if (consumesRef)
                    remaining -= firstOpLength;
            }

            // CIGAR straddles clip
            else {
                assert(firstOpLength > remaining);
                firstOp.Length(firstOpLength - remaining);
                if (consumesQuery)
                    queryPosRemovedFront += remaining;
                remaining = 0;
            }
        }
    }

    // -------------------------
    // clip trailing CIGAR ops
    // -------------------------

    remaining = origTEnd - newTEnd;
    while (remaining > 0 && !cigar.empty()) {
        CigarOperation& lastOp = cigar.back();
        const size_t lastOpLength = lastOp.Length();
        const bool consumesQuery = internal::ConsumesQuery(lastOp.Type());
        const bool consumesRef = internal::ConsumesReference(lastOp.Type());

        if (!consumesRef) {

            // e.g. softclip - just pop it completely
            cigar.pop_back();
            if (consumesQuery)
                queryPosRemovedBack += lastOpLength;

        } else {
            assert(consumesRef);

            // CIGAR ends at or after clip
            if (lastOpLength <= remaining) {
                cigar.pop_back();
                if (consumesQuery)
                    queryPosRemovedBack += lastOpLength;
                if (consumesRef)
                    remaining -= lastOpLength;
            }

            // CIGAR straddles clip
            else {
                assert(lastOpLength > remaining);
                lastOp.Length(lastOpLength - remaining);
                if (consumesQuery)
                    queryPosRemovedBack += remaining;
                remaining = 0;
            }
        }
    }

    // update CIGAR and position
    impl_.CigarData(cigar);
    impl_.Position(newTStart);

    // clip SEQ, QUAL, tags
    const Position qStart = origQStart + queryPosRemovedFront;
    const Position qEnd   = origQEnd   - queryPosRemovedBack;
    const size_t clipFrom = queryPosRemovedFront;
    const size_t clipLength = qEnd - qStart;
    ClipFields(clipFrom, clipLength);

    // update query start/end
    internal::CreateOrEdit(internal::tagName_queryStart, qStart, &impl_);
    internal::CreateOrEdit(internal::tagName_queryEnd,   qEnd,   &impl_);
//    UpdateName();

    // reset any cached aligned start/end
    ResetCachedPositions();
    return *this;
}

BamRecord& BamRecord::ClipToReferenceReverse(const PacBio::BAM::Position start,
                                             const PacBio::BAM::Position end)
{
    assert(IsMapped());
    assert(AlignedStrand() == Strand::REVERSE);

    // cache original coords
    const Position origQStart = QueryStart();
    const Position origQEnd   = QueryEnd();
    const Position origTStart = ReferenceStart();
    const Position origTEnd   = ReferenceEnd();

    // skip if already within requested clip range
    if (start <= origTStart && end >= origTEnd)
        return *this;
    assert(AlignedStart() >= origQStart);
    assert(AlignedEnd()   <= origQEnd);

    const Position newTStart = std::max(origTStart, start);
    const Position newTEnd   = std::min(origTEnd, end);

    Cigar cigar = std::move(impl_.CigarData());

    size_t queryPosRemovedFront = 0;
    size_t queryPosRemovedBack  = 0;

    // update CIGAR - clip front ops, then clip back ops
    size_t remaining = newTStart - origTStart;
    while (remaining > 0 && !cigar.empty()) {
        CigarOperation& firstOp = cigar.front();
        const CigarOperationType firstOpType = firstOp.Type();
        const size_t firstOpLength = firstOp.Length();
        const bool consumesQuery = internal::ConsumesQuery(firstOpType);
        const bool consumesRef = internal::ConsumesReference(firstOpType);

        if (!consumesRef) {

            // e.g. softclip - just pop it completely
            cigar.erase(cigar.begin());
            if (consumesQuery)
                queryPosRemovedBack += firstOpLength;

        } else {
            assert(consumesRef);

            // CIGAR ends at or before clip
            if (firstOpLength <= remaining) {
                cigar.erase(cigar.begin());
                if (consumesQuery)
                    queryPosRemovedBack += firstOpLength;
                if (consumesRef)
                    remaining -= firstOpLength;
            }

            // CIGAR straddles clip
            else {
                assert(firstOpLength > remaining);
                firstOp.Length(firstOpLength - remaining);
                if (consumesQuery)
                    queryPosRemovedBack += remaining;
                remaining = 0;
            }
        }
    }

    remaining = origTEnd - newTEnd;
    while (remaining > 0 && !cigar.empty()) {
        CigarOperation& lastOp = cigar.back();
        const CigarOperationType lastOpType = lastOp.Type();
        const size_t lastOpLength = lastOp.Length();
        const bool consumesQuery = internal::ConsumesQuery(lastOpType);
        const bool consumesRef = internal::ConsumesReference(lastOpType);

        if (!consumesRef) {

            // e.g. softclip - just pop it completely
            cigar.pop_back();
            if (consumesQuery)
                queryPosRemovedFront += lastOpLength;

        } else {
            assert(consumesRef);

            // CIGAR ends at or before clip
            if (lastOpLength <= remaining) {
                cigar.pop_back();
                if (consumesQuery)
                    queryPosRemovedFront += lastOpLength;
                if (consumesRef)
                    remaining -= lastOpLength;
            }

            // CIGAR straddles clip
            else {
                assert(lastOpLength > remaining);
                lastOp.Length(lastOpLength - remaining);
                if (consumesQuery)
                    queryPosRemovedFront += remaining;
                remaining = 0;
            }
        }
    }
    impl_.CigarData(cigar);

    // update aligned reference position
    impl_.Position(newTStart);

    // clip SEQ, QUAL, tags
    const Position qStart = origQStart + queryPosRemovedFront;
    const Position qEnd   = origQEnd   - queryPosRemovedBack;
    const size_t clipFrom = queryPosRemovedFront;
    const size_t clipLength = qEnd - qStart;
    ClipFields(clipFrom, clipLength);

    // update query start/end
    internal::CreateOrEdit(internal::tagName_queryStart, qStart, &impl_);
    internal::CreateOrEdit(internal::tagName_queryEnd,   qEnd,   &impl_);
//    UpdateName();

    // reset any cached aligned start/end
    ResetCachedPositions();
    return *this;
}

QualityValues BamRecord::DeletionQV(Orientation orientation,
                                    bool aligned,
                                    bool exciseSoftClips) const
{
    return FetchQualities(internal::tagName_deletionQV,
                          orientation,
                          aligned,
                          exciseSoftClips);
}

BamRecord& BamRecord::DeletionQV(const QualityValues& deletionQVs)
{
    internal::CreateOrEdit(internal::tagName_deletionQV, deletionQVs.Fastq(), &impl_);
    return *this;
}


string BamRecord::DeletionTag(Orientation orientation,
                              bool aligned,
                              bool exciseSoftClips) const
{
    return FetchBases(internal::tagName_deletionTag,
                      orientation,
                      aligned,
                      exciseSoftClips);
}

BamRecord& BamRecord::DeletionTag(const std::string& tags)
{
    internal::CreateOrEdit(internal::tagName_deletionTag, tags, &impl_);
    return *this;
}

std::vector<uint16_t>
BamRecord::EncodePhotons(const std::vector<float>& data)
{
    std::vector<uint16_t> encoded;
    encoded.reserve(data.size());
    for (const auto& d : data)
        encoded.emplace_back(d * photonFactor);
    return encoded;
}

string BamRecord::FetchBasesRaw(const string& tagName) const
{
    const Tag& seqTag = impl_.TagValue(tagName);
    return seqTag.ToString();
}

string BamRecord::FetchBases(const string& tagName,
                             const Orientation orientation) const
{ return FetchBases(tagName, orientation, false, false); }

string BamRecord::FetchBases(const string& tagName,
                             const Orientation orientation,
                             const bool aligned,
                             const bool exciseSoftClips) const
{
    // requested data info
    const bool isBamSeq = (tagName == internal::tagName_SEQ);
    const bool isPulse = (tagName == internal::tagName_pulse_call);

    // fetch raw
    string bases;
    Orientation current;
    if (isBamSeq) { // SEQ stored in genomic orientation
        bases = std::move(impl_.Sequence());
        current = Orientation::GENOMIC;
    } else { // all tags stored in native orientation
        bases = std::move(FetchBasesRaw(tagName));
        current = Orientation::NATIVE;
    }

    // if we need to touch CIGAR
    if (aligned || exciseSoftClips) {

        // force into genomic orientation
        internal::OrientBasesAsRequested(&bases,
                                         current,
                                         Orientation::GENOMIC,
                                         impl_.IsReverseStrand(),
                                         isPulse);
        current = Orientation::GENOMIC;

        // clip & gapify as requested
        internal::ClipAndGapifyBases(impl_,
                                     aligned,
                                     exciseSoftClips,
                                     &bases);
    }

    // return in the orientation requested
    internal::OrientBasesAsRequested(&bases,
                                     current,
                                     orientation,
                                     impl_.IsReverseStrand(),
                                     isPulse);
    return bases;
}

Frames BamRecord::FetchFramesRaw(const string& tagName) const
{
    Frames frames;
    const Tag& frameTag = impl_.TagValue(tagName);
    if (frameTag.IsNull())
        return frames;  // throw ?

    // lossy frame codes
    if (frameTag.IsUInt8Array()) {
        const vector<uint8_t> codes = std::move(frameTag.ToUInt8Array());
        frames = std::move(Frames::Decode(codes));
    }

    // lossless frame data
    else {
        assert(frameTag.IsUInt16Array());
        const vector<uint16_t> losslessFrames = std::move(frameTag.ToUInt16Array());
        frames.Data(std::move(losslessFrames));
    }

    return frames;
}

Frames BamRecord::FetchFrames(const string& tagName,
                              const Orientation orientation) const
{ return FetchFrames(tagName, orientation, false, false); }

Frames BamRecord::FetchFrames(const string& tagName,
                              const Orientation orientation,
                              const bool aligned,
                              const bool exciseSoftClips) const
{
    // fetch raw
    Frames frames = FetchFramesRaw(tagName);
    Orientation current = Orientation::NATIVE;

    if (aligned || exciseSoftClips) {

        // force into genomic orientation
        internal::OrientTagDataAsRequested(&frames,
                                           current,
                                           Orientation::GENOMIC,
                                           impl_.IsReverseStrand());
        current = Orientation::GENOMIC;

        // clip & gapify as requested
        internal::ClipAndGapifyFrames(impl_,
                                      aligned,
                                      exciseSoftClips,
                                      &frames);
    }

    // return in the orientation requested
    internal::OrientTagDataAsRequested(&frames,
                                       current,
                                       orientation,
                                       impl_.IsReverseStrand());
    return frames;
}

vector<float> BamRecord::FetchPhotons(const string& tagName,
                                      const Orientation orientation) const
{
    // fetch tag data
    const Tag& frameTag = impl_.TagValue(tagName);
    if (frameTag.IsNull())
        return vector<float>();
    if(!frameTag.IsUInt16Array())
        throw std::runtime_error("Photons are not a uint16_t array, tag " + tagName);
    vector<uint16_t> data = std::move(frameTag.ToUInt16Array());

    // put in requested orientation
    internal::OrientTagDataAsRequested(&data,
                                       Orientation::NATIVE,     // current
                                       orientation,             // requested
                                       impl_.IsReverseStrand());

    // store & return result
    vector<float> photons;
    photons.reserve(data.size());
    for (const auto& d : data)
        photons.emplace_back(d / photonFactor);
    return photons;
}

QualityValues BamRecord::FetchQualitiesRaw(const string& tagName) const
{
    const Tag& qvsTag = impl_.TagValue(tagName);
    return QualityValues::FromFastq(qvsTag.ToString());
}

QualityValues BamRecord::FetchQualities(const string& tagName,
                                        const Orientation orientation) const
{ return FetchQualities(tagName, orientation, false, false); }

QualityValues BamRecord::FetchQualities(const string& tagName,
                                        const Orientation orientation,
                                        const bool aligned,
                                        const bool exciseSoftClips) const
{
    // requested data info
    const bool isBamQual = (tagName == internal::tagName_QUAL);

    // fetch raw
    QualityValues quals;
    Orientation current;
    if (isBamQual) { // QUAL stored in genomic orientation
        quals = std::move(impl_.Qualities());
        current = Orientation::GENOMIC;
    } else {        // all tags stored in native orientation
        quals = std::move(FetchQualitiesRaw(tagName));
        current = Orientation::NATIVE;
    }

    // if we need to touch CIGAR
    if (aligned || exciseSoftClips) {

        // force into genomic orientation
        internal::OrientTagDataAsRequested(&quals,
                                          current,
                                          Orientation::GENOMIC,
                                          impl_.IsReverseStrand());
        current = Orientation::GENOMIC;

        // clip & gapify as requested
        internal::ClipAndGapifyQualities(impl_,
                                         aligned,
                                         exciseSoftClips,
                                         &quals);
    }

    // return in the orientation requested
    internal::OrientTagDataAsRequested(&quals,
                                       current,
                                       orientation,
                                       impl_.IsReverseStrand());
    return quals;
}

string BamRecord::FullName(void) const
{ return impl_.Name(); }

bool BamRecord::HasAltLabelQV(void) const
{ return impl_.HasTag(internal::tagName_alternative_labelQV); }

bool BamRecord::HasAltLabelTag(void) const
{ return impl_.HasTag(internal::tagName_alternative_labelTag); }

bool BamRecord::HasBarcodes(void) const
{ return impl_.HasTag(internal::tagName_barcodes); }

bool BamRecord::HasBarcodeQuality(void) const
{ return impl_.HasTag(internal::tagName_barcode_quality); }

bool BamRecord::HasLabelQV(void) const
{ return impl_.HasTag(internal::tagName_labelQV); }

bool BamRecord::HasDeletionQV(void) const
{ return impl_.HasTag(internal::tagName_deletionQV); }

bool BamRecord::HasDeletionTag(void) const
{ return impl_.HasTag(internal::tagName_deletionTag); }

bool BamRecord::HasHoleNumber(void) const
{ return impl_.HasTag(internal::tagName_holeNumber)
          && !impl_.TagValue(internal::tagName_holeNumber).IsNull();
}

bool BamRecord::HasInsertionQV(void) const
{ return impl_.HasTag(internal::tagName_insertionQV); }

bool BamRecord::HasNumPasses(void) const
{ return impl_.HasTag(internal::tagName_numPasses); }

bool BamRecord::HasPreBaseFrames(void) const
{ return HasIPD(); }

bool BamRecord::HasIPD(void) const
{ return impl_.HasTag(internal::tagName_ipd); }

bool BamRecord::HasLocalContextFlags(void) const
{ return impl_.HasTag(internal::tagName_contextFlags); }

bool BamRecord::HasMergeQV(void) const
{ return impl_.HasTag(internal::tagName_mergeQV); }

bool BamRecord::HasPulseMergeQV(void) const
{ return impl_.HasTag(internal::tagName_pulseMergeQV); }

bool BamRecord::HasPkmean(void) const
{ return impl_.HasTag(internal::tagName_pkmean); }

bool BamRecord::HasPkmid(void) const
{ return impl_.HasTag(internal::tagName_pkmid); }

bool BamRecord::HasPkmean2(void) const
{ return impl_.HasTag(internal::tagName_pkmean2); }

bool BamRecord::HasPkmid2(void) const
{ return impl_.HasTag(internal::tagName_pkmid2); }

bool BamRecord::HasPrePulseFrames(void) const
{ return impl_.HasTag(internal::tagName_pre_pulse_frames); }

bool BamRecord::HasPulseCall(void) const
{ return impl_.HasTag(internal::tagName_pulse_call)
          && !impl_.TagValue(internal::tagName_pulse_call).IsNull();
}

bool BamRecord::HasPulseCallWidth(void) const
{ return impl_.HasTag(internal::tagName_pulse_call_width); }

bool BamRecord::HasPulseWidth(void) const
{ return impl_.HasTag(internal::tagName_pulseWidth); }

bool BamRecord::HasQueryEnd(void) const
{ return impl_.HasTag(internal::tagName_queryEnd); }

bool BamRecord::HasQueryStart(void) const
{ return impl_.HasTag(internal::tagName_queryStart); }

bool BamRecord::HasReadAccuracy(void) const
{ return impl_.HasTag(internal::tagName_readAccuracy)
          && !impl_.TagValue(internal::tagName_readAccuracy).IsNull();
}

bool BamRecord::HasScrapRegionType(void) const
{ return impl_.HasTag(internal::tagName_scrap_region_type)
          && !impl_.TagValue(internal::tagName_scrap_region_type).IsNull();
}

bool BamRecord::HasScrapZmwType(void) const
{ return impl_.HasTag(internal::tagName_scrap_zmw_type)
          && !impl_.TagValue(internal::tagName_scrap_zmw_type).IsNull();
}

bool BamRecord::HasStartFrame(void) const
{ return impl_.HasTag(internal::tagName_startFrame); }

bool BamRecord::HasSignalToNoise(void) const
{ return impl_.HasTag(internal::tagName_snr); }

bool BamRecord::HasSubstitutionQV(void) const
{ return impl_.HasTag(internal::tagName_substitutionQV); }

bool BamRecord::HasSubstitutionTag(void) const
{ return impl_.HasTag(internal::tagName_substitutionTag); }

BamHeader BamRecord::Header(void) const
{ return header_; }

int32_t BamRecord::HoleNumber(void) const
{
    const Tag& holeNumber = impl_.TagValue(internal::tagName_holeNumber);
    if (!holeNumber.IsNull())
        return holeNumber.ToInt32();

    // missing zm tag - try to pull from name
    return internal::HoleNumberFromName(FullName());
}

BamRecord& BamRecord::HoleNumber(const int32_t holeNumber)
{
    internal::CreateOrEdit(internal::tagName_holeNumber,
                           holeNumber,
                           &impl_);
    return *this;
}

BamRecordImpl& BamRecord::Impl(void)
{ return impl_; }

const BamRecordImpl& BamRecord::Impl(void) const
{ return impl_; }

QualityValues BamRecord::InsertionQV(Orientation orientation,
                                     bool aligned,
                                     bool exciseSoftClips) const
{
    return FetchQualities(internal::tagName_insertionQV,
                          orientation,
                          aligned,
                          exciseSoftClips);
}

BamRecord& BamRecord::InsertionQV(const QualityValues& insertionQVs)
{
    internal::CreateOrEdit(internal::tagName_insertionQV, insertionQVs.Fastq(), &impl_);
    return *this;
}

Frames BamRecord::IPD(Orientation orientation,
                      bool aligned,
                      bool exciseSoftClips) const
{
    return FetchFrames(internal::tagName_ipd,
                       orientation,
                       aligned,
                       exciseSoftClips);
}

BamRecord& BamRecord::IPD(const Frames& frames,
                          const FrameEncodingType encoding)
{
    if (encoding == FrameEncodingType::LOSSY)
        internal::CreateOrEdit(internal::tagName_ipd, frames.Encode(), &impl_);
    else
        internal::CreateOrEdit(internal::tagName_ipd, frames.Data(), &impl_);
    return *this;
}

size_t BamRecord::NumDeletedBases(void) const
{
    auto tEnd = ReferenceEnd();
    auto tStart = ReferenceStart();
    auto numMatchesAndMismatches = NumMatchesAndMismatches();
    auto nM = numMatchesAndMismatches.first;
    auto nMM = numMatchesAndMismatches.second;
    return (tEnd - tStart - nM - nMM);
}

size_t BamRecord::NumInsertedBases(void) const
{
    auto aEnd = AlignedEnd();
    auto aStart = AlignedStart();
    auto numMatchesAndMismatches = NumMatchesAndMismatches();
    auto nM = numMatchesAndMismatches.first;
    auto nMM = numMatchesAndMismatches.second;
    return (aEnd - aStart - nM - nMM);
}

size_t BamRecord::NumMatches(void) const
{ return NumMatchesAndMismatches().first; }

pair<size_t, size_t> BamRecord::NumMatchesAndMismatches(void) const
{
    pair<size_t, size_t> result = make_pair(0,0);
    PBBAM_SHARED_PTR<bam1_t> b = internal::BamRecordMemory::GetRawData(this);
    uint32_t* cigarData = bam_get_cigar(b.get());
    for (uint32_t i = 0; i < b->core.n_cigar; ++i) {
        const CigarOperationType type = static_cast<CigarOperationType>(bam_cigar_op(cigarData[i]));
        if (type == CigarOperationType::SEQUENCE_MATCH)
            result.first += bam_cigar_oplen(cigarData[i]);
        else if (type == CigarOperationType::SEQUENCE_MISMATCH)
            result.second += bam_cigar_oplen(cigarData[i]);
    }
    return result;
}

size_t BamRecord::NumMismatches(void) const
{ return NumMatchesAndMismatches().second; }

Frames BamRecord::PreBaseFrames(Orientation orientation, 
                                bool aligned,
                                bool exciseSoftClips) const
{ return IPD(orientation, aligned, exciseSoftClips); }

BamRecord& BamRecord::PreBaseFrames(const Frames& frames,
                                    const FrameEncodingType encoding)
{ return IPD(frames, encoding); }

Frames BamRecord::IPDRaw(Orientation orientation) const
{
    const auto tagName = internal::tagName_ipd;

    Frames frames;
    const Tag& frameTag = impl_.TagValue(tagName);
    if (frameTag.IsNull())
        return frames;

    // lossy frame codes
    if (frameTag.IsUInt8Array()) {
        const vector<uint8_t> codes = std::move(frameTag.ToUInt8Array());
        const vector<uint16_t> codes16(codes.begin(), codes.end());
        frames.Data(std::move(codes16));
    }

    // lossless frame data
    else {
        assert(frameTag.IsUInt16Array());
        const vector<uint16_t> losslessFrames = std::move(frameTag.ToUInt16Array());
        frames.Data(std::move(losslessFrames));
    }

    // return in requested orientation
    internal::OrientTagDataAsRequested(&frames,
                                       Orientation::NATIVE,     // current
                                       orientation,             // requested
                                       impl_.IsReverseStrand());
    return frames;
}

Frames BamRecord::PulseWidthRaw(Orientation orientation) const
{
    const auto tagName = internal::tagName_pulseWidth;

    Frames frames;
    const Tag& frameTag = impl_.TagValue(tagName);
    if (frameTag.IsNull())
        return frames;

    // lossy frame codes
    if (frameTag.IsUInt8Array()) {
        const vector<uint8_t> codes = std::move(frameTag.ToUInt8Array());
        const vector<uint16_t> codes16(codes.begin(), codes.end());
        frames.Data(std::move(codes16));
    }

    // lossless frame data
    else {
        assert(frameTag.IsUInt16Array());
        const vector<uint16_t> losslessFrames = std::move(frameTag.ToUInt16Array());
        frames.Data(std::move(losslessFrames));
    }

    // return in requested orientation
    internal::OrientTagDataAsRequested(&frames,
                                       Orientation::NATIVE,  // current
                                       orientation,          // requested
                                       impl_.IsReverseStrand());
    return frames;
}

bool BamRecord::IsMapped(void) const
{ return impl_.IsMapped(); }

QualityValues BamRecord::LabelQV(Orientation orientation) const
{ return FetchQualities(internal::tagName_labelQV, orientation); }

BamRecord& BamRecord::LabelQV(const QualityValues& labelQVs)
{
    internal::CreateOrEdit(internal::tagName_labelQV, labelQVs.Fastq(), &impl_);
    return *this;
}

LocalContextFlags BamRecord::LocalContextFlags(void) const
{
    const Tag& cxTag = impl_.TagValue(internal::tagName_contextFlags);
    return static_cast<PacBio::BAM::LocalContextFlags>(cxTag.ToUInt8());
}

BamRecord& BamRecord::LocalContextFlags(const PacBio::BAM::LocalContextFlags flags)
{
    internal::CreateOrEdit(internal::tagName_contextFlags,
                           static_cast<uint8_t>(flags),
                           &impl_);
    return *this;
}

BamRecord& BamRecord::Map(const int32_t referenceId,
                          const Position refStart,
                          const Strand strand,
                          const Cigar& cigar,
                          const uint8_t mappingQuality)
{
    impl_.Position(refStart);
    impl_.ReferenceId(referenceId);
    impl_.CigarData(cigar);
    impl_.MapQuality(mappingQuality);
    impl_.SetMapped(true);

    if (strand == Strand::FORWARD)
        impl_.SetReverseStrand(false);

    else {
        assert(strand == Strand::REVERSE);
        impl_.SetReverseStrand(true);

        // switch seq & qual
        string sequence  = impl_.Sequence();
        QualityValues qualities = impl_.Qualities();

        internal::ReverseComplement(sequence);
        internal::Reverse(qualities);

        impl_.SetSequenceAndQualities(sequence, qualities.Fastq());
    }

    // reset any cached aligned start/end
    alignedStart_ = PacBio::BAM::UnmappedPosition;
    alignedEnd_ = PacBio::BAM::UnmappedPosition;

    return *this;
}

uint8_t BamRecord::MapQuality(void) const
{ return impl_.MapQuality(); }

QualityValues BamRecord::MergeQV(Orientation orientation,
                                 bool aligned,
                                 bool exciseSoftClips) const
{
    return FetchQualities(internal::tagName_mergeQV,
                          orientation,
                          aligned,
                          exciseSoftClips);
}

BamRecord& BamRecord::MergeQV(const QualityValues& mergeQVs)
{
    internal::CreateOrEdit(internal::tagName_mergeQV, mergeQVs.Fastq(), &impl_);
    return *this;
}

QualityValues BamRecord::PulseMergeQV(Orientation orientation) const
{ return FetchQualities(internal::tagName_pulseMergeQV, orientation); }

BamRecord& BamRecord::PulseMergeQV(const QualityValues& mergeQVs)
{
    internal::CreateOrEdit(internal::tagName_pulseMergeQV, mergeQVs.Fastq(), &impl_);
    return *this;
}

string BamRecord::MovieName(void) const
{ return ReadGroup().MovieName(); }

int32_t BamRecord::NumPasses(void) const
{
    const Tag& numPasses = impl_.TagValue(internal::tagName_numPasses);
    return numPasses.ToInt32();
}

BamRecord& BamRecord::NumPasses(const int32_t numPasses)
{
    internal::CreateOrEdit(internal::tagName_numPasses, numPasses, &impl_);
    return *this;
}

std::vector<float> BamRecord::Pkmean(Orientation orientation) const
{ return FetchPhotons(internal::tagName_pkmean, orientation); }

BamRecord& BamRecord::Pkmean(const std::vector<float>& photons)
{
    Pkmean(EncodePhotons(photons));
    return *this;
}

BamRecord& BamRecord::Pkmean(const std::vector<uint16_t>& encodedPhotons)
{
    internal::CreateOrEdit(internal::tagName_pkmean, encodedPhotons, &impl_);
    return *this;
}

std::vector<float> BamRecord::Pkmid(Orientation orientation) const
{ return FetchPhotons(internal::tagName_pkmid, orientation); }

BamRecord& BamRecord::Pkmid(const std::vector<float>& photons)
{
    Pkmid(EncodePhotons(photons));
    return *this;
}

BamRecord& BamRecord::Pkmid(const std::vector<uint16_t>& encodedPhotons)
{
    internal::CreateOrEdit(internal::tagName_pkmid, encodedPhotons, &impl_);
    return *this;
}

std::vector<float> BamRecord::Pkmean2(Orientation orientation) const
{ return FetchPhotons(internal::tagName_pkmean2, orientation); }

BamRecord& BamRecord::Pkmean2(const std::vector<float>& photons)
{
    Pkmean2(EncodePhotons(photons));
    return *this;
}

BamRecord& BamRecord::Pkmean2(const std::vector<uint16_t>& encodedPhotons)
{
    internal::CreateOrEdit(internal::tagName_pkmean2, encodedPhotons, &impl_);
    return *this;
}

std::vector<float> BamRecord::Pkmid2(Orientation orientation) const
{ return FetchPhotons(internal::tagName_pkmid2, orientation); }

BamRecord& BamRecord::Pkmid2(const std::vector<float>& photons)
{
    Pkmid2(EncodePhotons(photons));
    return *this;
}

BamRecord& BamRecord::Pkmid2(const std::vector<uint16_t>& encodedPhotons)
{
    internal::CreateOrEdit(internal::tagName_pkmid2, encodedPhotons, &impl_);
    return *this;
}

Frames BamRecord::PrePulseFrames(Orientation orientation) const
{ return FetchFrames(internal::tagName_pre_pulse_frames, orientation); }

BamRecord& BamRecord::PrePulseFrames(const Frames& frames,
                                     const FrameEncodingType encoding)
{
    if (encoding == FrameEncodingType::LOSSY)
        internal::CreateOrEdit(internal::tagName_pre_pulse_frames, frames.Encode(), &impl_);
    else
        internal::CreateOrEdit(internal::tagName_pre_pulse_frames, frames.Data(), &impl_);
    return *this;
}

std::string BamRecord::PulseCall(Orientation orientation) const
{ return FetchBases(internal::tagName_pulse_call, orientation); }

BamRecord& BamRecord::PulseCall(const std::string& tags)
{
    internal::CreateOrEdit(internal::tagName_pulse_call, tags, &impl_);
    return *this;
}

Frames BamRecord::PulseCallWidth(Orientation orientation) const
{ return FetchFrames(internal::tagName_pulse_call_width, orientation); }

BamRecord& BamRecord::PulseCallWidth(const Frames& frames,
                                     const FrameEncodingType encoding)
{
    if (encoding == FrameEncodingType::LOSSY)
        internal::CreateOrEdit(internal::tagName_pulse_call_width, frames.Encode(), &impl_);
    else
        internal::CreateOrEdit(internal::tagName_pulse_call_width, frames.Data(), &impl_);
    return *this;
}

Frames BamRecord::PulseWidth(Orientation orientation,
                             bool aligned,
                             bool exciseSoftClips) const
{
    return FetchFrames(internal::tagName_pulseWidth,
                       orientation,
                       aligned,
                       exciseSoftClips);
}

BamRecord& BamRecord::PulseWidth(const Frames& frames,
                                 const FrameEncodingType encoding)
{
    if (encoding == FrameEncodingType::LOSSY)
        internal::CreateOrEdit(internal::tagName_pulseWidth, frames.Encode(), &impl_);
    else
        internal::CreateOrEdit(internal::tagName_pulseWidth, frames.Data(), &impl_);
    return *this;
}

QualityValues BamRecord::Qualities(Orientation orientation,
                                   bool aligned,
                                   bool exciseSoftClips) const
{
    return FetchQualities("QUAL",
                          orientation,
                          aligned,
                          exciseSoftClips);
}

Position BamRecord::QueryEnd(void) const
{
    // try 'qe' tag
    const Tag& qe = impl_.TagValue(internal::tagName_queryEnd);
    if (!qe.IsNull())
        return qe.ToInt32();

    // tag missing, need to check movie name (fallback for non-PB BAMs, but ignore for CCS reads)
    RecordType type;
    try {
        type = Type();
    } catch (std::exception&) {
        return Position(0);
    }
    if (type == RecordType::CCS)
        throw std::runtime_error("no query end for CCS read type");

    // PacBio BAM, non-CCS
    try {
        return internal::QueryEndFromName(FullName());
    } catch (std::exception&) {
        // return fallback position
        return Position(0);
    }
}

BamRecord& BamRecord::QueryEnd(const Position pos)
{
   internal::CreateOrEdit(internal::tagName_queryEnd,
                          (int32_t)pos,
                          &impl_);
   UpdateName();
   return *this;
}

Position BamRecord::QueryStart(void) const
{
    // try 'qs' tag
    const Tag& qs = impl_.TagValue(internal::tagName_queryStart);
    if (!qs.IsNull())
        return qs.ToInt32();

    // tag missing, need to check movie name (fallback for non-PB BAMs, but ignore for CCS reads)
    RecordType type;
    try {
        type = Type();
    } catch (std::exception&) {
        return Position(0);
    }
    if (type == RecordType::CCS)
        throw std::runtime_error("no query start for CCS read type");

    // PacBio BAM, non-CCS
    try {
        return internal::QueryStartFromName(FullName());
    } catch (std::exception&) {
        // return fallback position
        return Position(0);
    }
}

BamRecord& BamRecord::QueryStart(const Position pos)
{
   internal::CreateOrEdit(internal::tagName_queryStart,
                          (int32_t)pos,
                          &impl_);
   UpdateName();
   return *this;
}

Accuracy BamRecord::ReadAccuracy(void) const
{
    const Tag& readAccuracy = impl_.TagValue(internal::tagName_readAccuracy);
    return Accuracy(readAccuracy.ToFloat());
}

BamRecord& BamRecord::ReadAccuracy(const Accuracy& accuracy)
{
    internal::CreateOrEdit(internal::tagName_readAccuracy,
                           static_cast<float>(accuracy),
                           &impl_);
    return *this;
}

ReadGroupInfo BamRecord::ReadGroup(void) const
{ return header_.ReadGroup(ReadGroupId()); }

BamRecord& BamRecord::ReadGroup(const ReadGroupInfo& rg)
{
   internal::CreateOrEdit(internal::tagName_readGroup, rg.Id(), &impl_);
   UpdateName();
   return *this;
}

string BamRecord::ReadGroupId(void) const
{
    const Tag& rgTag = impl_.TagValue(internal::tagName_readGroup);
    if (rgTag.IsNull())
        return string();
    return rgTag.ToString();
}

BamRecord& BamRecord::ReadGroupId(const std::string& id)
{
   internal::CreateOrEdit(internal::tagName_readGroup,
                          id,
                          &impl_);
   UpdateName();
   return *this;
}

int32_t BamRecord::ReadGroupNumericId(void) const
{ return ReadGroupInfo::IdToInt(ReadGroupId()); }

Position BamRecord::ReferenceEnd(void) const
{
    if (!impl_.IsMapped())
        return PacBio::BAM::UnmappedPosition;
    PBBAM_SHARED_PTR<bam1_t> htsData = internal::BamRecordMemory::GetRawData(impl_);
    if (!htsData)
        return PacBio::BAM::UnmappedPosition;
    return bam_endpos(htsData.get());
}

int32_t BamRecord::ReferenceId(void) const
{ return impl_.ReferenceId(); }

std::string BamRecord::ReferenceName(void) const
{
    if (IsMapped())
        return Header().SequenceName(ReferenceId());
    else
        throw std::runtime_error("unmapped record has no associated reference name");
}

Position BamRecord::ReferenceStart(void) const
{ return impl_.Position(); }

void BamRecord::ResetCachedPositions(void) const
{
    alignedEnd_   = PacBio::BAM::UnmappedPosition;
    alignedStart_ = PacBio::BAM::UnmappedPosition;
}

void BamRecord::ResetCachedPositions(void)
{
    alignedEnd_   = PacBio::BAM::UnmappedPosition;
    alignedStart_ = PacBio::BAM::UnmappedPosition;
}

VirtualRegionType BamRecord::ScrapRegionType(void) const
{
    const Tag& srTag = impl_.TagValue(internal::tagName_scrap_region_type);
    return VirtualRegionTypeMap::ParseChar[srTag.ToUInt8()];
}

BamRecord& BamRecord::ScrapRegionType(const VirtualRegionType type)
{
    internal::CreateOrEdit(internal::tagName_scrap_region_type,
                           static_cast<uint8_t>(type), &impl_);
    return *this;
}

BamRecord& BamRecord::ScrapRegionType(const char type)
{
    internal::CreateOrEdit(internal::tagName_scrap_region_type, type, &impl_);
    return *this;
}

ZmwType BamRecord::ScrapZmwType(void) const
{
    const Tag& szTag = impl_.TagValue(internal::tagName_scrap_zmw_type);
    return ZmwTypeMap::ParseChar[szTag.ToUInt8()];
}

BamRecord& BamRecord::ScrapZmwType(const ZmwType type)
{
    internal::CreateOrEdit(internal::tagName_scrap_zmw_type,
                           static_cast<uint8_t>(type), &impl_);
    return *this;
}

BamRecord& BamRecord::ScrapZmwType(const char type)
{
    internal::CreateOrEdit(internal::tagName_scrap_zmw_type, type, &impl_);
    return *this;
}

std::string BamRecord::Sequence(const Orientation orientation,
                                bool aligned,
                                bool exciseSoftClips) const
{
    return FetchBases("SEQ",
                      orientation,
                      aligned,
                      exciseSoftClips);
}

vector<float> BamRecord::SignalToNoise(void) const
{
    const Tag& snTag = impl_.TagValue(internal::tagName_snr);
    return snTag.ToFloatArray();
}

BamRecord& BamRecord::SignalToNoise(const vector<float>& snr)
{
    internal::CreateOrEdit(internal::tagName_snr, snr, &impl_);
    return *this;
}

std::vector<uint32_t> BamRecord::StartFrame(Orientation orientation) const
{
    const Tag& sfTag = impl_.TagValue(internal::tagName_startFrame);
    return sfTag.ToUInt32Array();
}

BamRecord& BamRecord::StartFrame(const std::vector<uint32_t>& startFrame)
{
    internal::CreateOrEdit(internal::tagName_startFrame, startFrame, &impl_);
    return *this;
}

QualityValues BamRecord::SubstitutionQV(Orientation orientation,
                                        bool aligned,
                                        bool exciseSoftClips) const
{
    return FetchQualities(internal::tagName_substitutionQV,
                          orientation,
                          aligned,
                          exciseSoftClips);
}

BamRecord& BamRecord::SubstitutionQV(const QualityValues& substitutionQVs)
{
    internal::CreateOrEdit(internal::tagName_substitutionQV, substitutionQVs.Fastq(), &impl_);
    return *this;
}

std::string BamRecord::SubstitutionTag(Orientation orientation,
                                        bool aligned,
                                        bool exciseSoftClips) const
{
    return FetchBases(internal::tagName_substitutionTag,
                      orientation,
                      aligned,
                      exciseSoftClips);
}

BamRecord& BamRecord::SubstitutionTag(const std::string& tags)
{
    internal::CreateOrEdit(internal::tagName_substitutionTag, tags, &impl_);
    return *this;
}

RecordType BamRecord::Type(void) const
{
    try {
        const string& typeName = ReadGroup().ReadType();
        return internal::NameToType(typeName);
    } catch (std::exception&) {

        // read group not found
        // peek at name to see if we're CCS
        if (FullName().find("ccs") != string::npos)
            return RecordType::CCS;

        // otherwise unknown
        else
            return RecordType::UNKNOWN;
    }
}

void BamRecord::UpdateName()
{
    std::string newName;
    newName.reserve(100);

    newName += MovieName();
    newName += "/";

    if (HasHoleNumber())
        newName += std::to_string(HoleNumber());
    else
        newName += "?";

    newName += "/";

    if (Type() == RecordType::CCS)
        newName += "ccs";
    else {
        if (HasQueryStart())
            newName += std::to_string(QueryStart());
        else
            newName += "?";

        newName += '_';

        if (HasQueryEnd())
            newName += std::to_string(QueryEnd());
        else
            newName += "?";
    }

    impl_.Name(newName);
}
