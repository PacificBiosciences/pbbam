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

#include "pbbam/BamRecord.h"
#include "AssertUtils.h"
#include "MemoryUtils.h"
#include "SequenceUtils.h"
#include <htslib/sam.h>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

// BAM record tag names
static const string tagName_readAccuracy    = "rq";
static const string tagName_holeNumber      = "zm";
static const string tagName_numPasses       = "np";
static const string tagName_deletionQV      = "dq";
static const string tagName_deletionTag     = "dt";
static const string tagName_insertionQV     = "iq";
static const string tagName_ipd             = "ip";
static const string tagName_mergeQV         = "mq";
static const string tagName_pulseWidth      = "pw";
static const string tagName_readGroup       = "RG";
static const string tagName_queryStart      = "qs";
static const string tagName_queryEnd        = "qe";
static const string tagName_substitutionQV  = "sq";
static const string tagName_substitutionTag = "st";

// faux (helper) tag names
static const string tagName_QUAL = "QUAL";
static const string tagName_SEQ  = "SEQ";

// record type names
static const string recordTypeName_Polymerase = "POLYMERASE";
static const string recordTypeName_HqRegion   = "HQREGION";
static const string recordTypeName_Subread    = "SUBREAD";
static const string recordTypeName_CCS        = "CCS";
static const string recordTypeName_Scrap      = "SCRAP";
static const string recordTypeName_Unknown    = "UNKNOWN";

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
int32_t AlignedEndOffset(const Cigar& cigar,
                         const int seqLength)
{
    int32_t endOffset = seqLength;

    if (!cigar.empty()) {
        Cigar::const_reverse_iterator cigarIter = cigar.crbegin();
        Cigar::const_reverse_iterator cigarEnd  = cigar.crend();
        for (; cigarIter != cigarEnd; ++cigarIter) {
            const CigarOperation& op = (*cigarIter);
            if (op.Type() == CigarOperationType::HARD_CLIP) {
                if (endOffset != 0 && endOffset != seqLength)
                    return -1;
            }
            else if (op.Type() == CigarOperationType::SOFT_CLIP)
                endOffset -= op.Length();
            else
                break;
        }
    }

    if (endOffset == 0)
        endOffset = seqLength;
    return endOffset;
}

static
int32_t AlignedStartOffset(const Cigar& cigar,
                           const int seqLength)
{
    int32_t startOffset = 0;

    if (!cigar.empty()) {
        Cigar::const_iterator cigarIter = cigar.cbegin();
        Cigar::const_iterator cigarEnd  = cigar.cend();
        for (; cigarIter != cigarEnd; ++cigarIter) {
            const CigarOperation& op = (*cigarIter);
            if (op.Type() == CigarOperationType::HARD_CLIP) {
                if (startOffset != 0 && startOffset != seqLength)
                    return -1;
            }
            else if (op.Type() == CigarOperationType::SOFT_CLIP)
                startOffset += op.Length();
            else
                break;
        }
    }
    return startOffset;
}

template<typename T>
T Clip(const T& input,
       const size_t pos,
       const size_t len)
{
    return T(input.cbegin() + pos,
             input.cbegin() + pos + len);
}


static
void MaybeClipAndGapifyBases(const BamRecordImpl& impl,
                             const bool aligned,
                             const bool exciseSoftClips,
                             string& seq)
{
    if (impl.IsMapped() && (aligned || exciseSoftClips)) {

        size_t seqIndex = 0;
        const Cigar& cigar = impl.CigarData();
        Cigar::const_iterator cigarIter = cigar.cbegin();
        Cigar::const_iterator cigarEnd  = cigar.cend();
        for (; cigarIter != cigarEnd; ++cigarIter) {
            const CigarOperation& op = (*cigarIter);
            const CigarOperationType& type = op.Type();

            // do nothing for hard clips
            if (type != CigarOperationType::HARD_CLIP) {
                const size_t opLength = op.Length();

                // maybe remove soft clips
                if (type == CigarOperationType::SOFT_CLIP && exciseSoftClips)
                    seq.erase(seqIndex, opLength);

                // for non-clipping operations
                else {

                    // maybe add gaps/padding
                    if (aligned) {
                        if (type == CigarOperationType::DELETION)
                            seq.insert(seqIndex, opLength, '-');
                        else if (type == CigarOperationType::PADDING)
                            seq.insert(seqIndex, opLength, '*');
                    }

                    // update index
                    seqIndex += opLength;
                }
            }
        }
    }
}

static
void MaybeClipAndGapifyFrames(const BamRecordImpl& impl,
                              const bool aligned,
                              const bool exciseSoftClips,
                              Frames& frames)
{
    if (impl.IsMapped() && (aligned || exciseSoftClips)) {

        vector<uint16_t> data = std::move(frames.Data()); // we're going to put it back
        size_t frameIndex = 0;
        const Cigar& cigar = impl.CigarData();
        Cigar::const_iterator cigarIter = cigar.cbegin();
        Cigar::const_iterator cigarEnd  = cigar.cend();
        for (; cigarIter != cigarEnd; ++cigarIter) {
            const CigarOperation& op = (*cigarIter);
            const CigarOperationType& type = op.Type();

            // do nothing for hard clips
            if (type != CigarOperationType::HARD_CLIP) {
                const size_t opLength = op.Length();

                // maybe remove soft clips
                if (type == CigarOperationType::SOFT_CLIP && exciseSoftClips)
                    data.erase(data.begin() + frameIndex, data.begin() + frameIndex + opLength);

                // for non-clipping operations
                else {

                    // maybe add gaps/padding
                    if (aligned) {
                        if (type == CigarOperationType::DELETION || type == CigarOperationType::PADDING)
                            data.insert(data.begin() + frameIndex, opLength, 0);
                    }

                    // update index
                    frameIndex += opLength;
                }
            }
        }
        frames.Data(data);
    }
}

static
void MaybeClipAndGapifyQualities(const BamRecordImpl& impl,
                                 const bool aligned,
                                 const bool exciseSoftClips,
                                 QualityValues& quals)
{
    if (impl.IsMapped() && (aligned || exciseSoftClips)) {

        size_t qualIndex = 0;
        const Cigar& cigar = impl.CigarData();
        Cigar::const_iterator cigarIter = cigar.cbegin();
        Cigar::const_iterator cigarEnd  = cigar.cend();
        for (; cigarIter != cigarEnd; ++cigarIter) {
            const CigarOperation& op = (*cigarIter);
            const CigarOperationType& type = op.Type();

            // do nothing for hard clips
            if (type != CigarOperationType::HARD_CLIP) {
                const size_t opLength = op.Length();

                // maybe remove soft clips
                if (type == CigarOperationType::SOFT_CLIP && exciseSoftClips)
                    quals.erase(quals.begin() + qualIndex, quals.begin() + qualIndex + opLength);

                // for non-clipping operations
                else {

                    // maybe add gaps/padding
                    if (aligned) {
                        if (type == CigarOperationType::DELETION || type == CigarOperationType::PADDING)
                            quals.insert(quals.begin() + qualIndex, opLength, QualityValue(0));
                    }

                    // update index
                    qualIndex += opLength;
                }
            }
        }
    }
}

static inline
void MaybeReverseFrames(const bool isReverseStrand,
                        const Orientation orientation,
                        Frames& frames)
{
    const bool shouldReverse = isReverseStrand && orientation == Orientation::GENOMIC;
    if (shouldReverse)
        std::reverse(frames.begin(), frames.end());
}

static inline
void MaybeReverseQuals(const bool isBamQual,
                       const bool isReverseStrand,
                       const Orientation orientation,
                       QualityValues& quals)
{
    const bool shouldReverse = (isBamQual ? isReverseStrand && orientation == Orientation::NATIVE
                                          : isReverseStrand && orientation == Orientation::GENOMIC);
    if (shouldReverse)
        std::reverse(quals.begin(), quals.end());
}

static inline
void MaybeReverseComplementSeq(const bool isBamSeq,
                               const bool isReverseStrand,
                               const Orientation orientation,
                               string& seq)
{
    const bool shouldReverse = (isBamSeq ? isReverseStrand && orientation == Orientation::NATIVE
                                         : isReverseStrand && orientation == Orientation::GENOMIC);
    if (shouldReverse)
        internal::ReverseComplement(seq);
}

//static
//string TypeToName(const RecordType& type)
//{
//    switch(type)
//    {
//        case RecordType::POLYMERASE : return recordTypeName_Polymerase;
//        case RecordType::HQREGION   : return recordTypeName_HqRegion;
//        case RecordType::SUBREAD    : return recordTypeName_Subread;
//        case RecordType::CCS        : return recordTypeName_CCS;
//        case RecordType::SCRAP      : return recordTypeName_Scrap;
//        case RecordType::UNKNOWN    : return recordTypeName_Unknown;
//        default:
//            PB_ASSERT_OR_RETURN_VALUE(false, string());
//    }
//}

static
RecordType NameToType(const string& name)
{
    if (name == recordTypeName_Subread)
        return RecordType::SUBREAD;
    if (name == recordTypeName_Polymerase)
        return RecordType::POLYMERASE;
    if (name == recordTypeName_HqRegion)
        return RecordType::HQREGION;
    if (name == recordTypeName_CCS)
        return RecordType::CCS;
    if (name == recordTypeName_Scrap)
        return RecordType::SCRAP;
    return RecordType::UNKNOWN;
}


} // namespace internal
} // namespace BAM
} // namespace PacBio

BamRecord::BamRecord(void)
    : header_(nullptr)
    , alignedStart_(PacBio::BAM::UnmappedPosition)
    , alignedEnd_(PacBio::BAM::UnmappedPosition)
{ }

BamRecord::BamRecord(const BamHeader::SharedPtr &header)
    : header_(header)
    , alignedStart_(PacBio::BAM::UnmappedPosition)
    , alignedEnd_(PacBio::BAM::UnmappedPosition)
{ }

BamRecord::BamRecord(const BamRecordImpl& impl)
    : impl_(impl)
    , header_(nullptr)
    , alignedStart_(PacBio::BAM::UnmappedPosition)
    , alignedEnd_(PacBio::BAM::UnmappedPosition)
{ }

BamRecord::BamRecord(BamRecordImpl&& impl)
    : impl_(std::move(impl))
    , header_(nullptr)
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

void BamRecord::CalculateAlignedPositions(void) const
{
    // reset
    alignedEnd_   = PacBio::BAM::UnmappedPosition;
    alignedStart_ = PacBio::BAM::UnmappedPosition;

    // skip if unmapped, or has no queryStart/End
    if (!impl_.IsMapped())
        return;
    const Position qStart = QueryStart();
    const Position qEnd   = QueryEnd();
    if (qStart == PacBio::BAM::UnmappedPosition || qEnd == PacBio::BAM::UnmappedPosition)
        return;

    // determine clipped end ranges
    const Cigar& cigar     = impl_.CigarData();
    const size_t seqLength = impl_.Sequence().size();
    const int32_t startOffset = internal::AlignedStartOffset(cigar, seqLength);
    const int32_t endOffset   = internal::AlignedEndOffset(cigar, seqLength);
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

Cigar BamRecord::CigarData(void) const
{ return impl_.CigarData(); }

BamRecord& BamRecord::Clip(const ClipType clipType,
                           const Position start,
                           const Position end)
{
    // skip if no clip requested
    if (clipType == ClipType::CLIP_NONE)
        return *this;
    const bool clipToQuery = (clipType == ClipType::CLIP_TO_QUERY);

    // cache original coords
    const Position origQStart = QueryStart();
    const Position origQEnd   = QueryEnd();
    const Position origAStart = AlignedStart();
    const Position origAEnd   = AlignedEnd();

    // only used on mapped records
    Position origTStart;
    Position origTEnd;
    bool isForwardStrand = (AlignedStrand() == Strand::FORWARD);

    // cache any add'l coords, skip out if clip not needed (or not possible)
    if (clipToQuery) {
        if (start <= origQStart && end >= origQEnd)
            return *this;
    } else {

        assert(clipType == ClipType::CLIP_TO_REFERENCE);
        if (!IsMapped())
            return *this;

        origTStart = ReferenceStart();
        origTEnd   = ReferenceEnd();
        if (start <= origTStart && end >= origTEnd)
            return *this;

        assert(origAStart >= origQStart);
        assert(origAEnd   <= origQEnd);
    }

    // determine new offsets into data
    size_t startOffset;
    size_t endOffset;

    if (clipToQuery) {
        startOffset = start - origQStart;
        endOffset   = origQEnd - end;
    } else {

        const size_t alignedStartOffset = (origAStart - origQStart);
        const size_t alignedEndOffset = (origQEnd - origAEnd);
        const size_t tStartDiff = start - origTStart;
        const size_t tEndDiff   = origTEnd - end;

        if (isForwardStrand) {
            startOffset = alignedStartOffset + tStartDiff;
            endOffset = alignedEndOffset + tEndDiff;
        } else {
            startOffset = alignedEndOffset + tStartDiff;
            endOffset = alignedStartOffset + tEndDiff;
        }
    }

    size_t queryPosRemovedFront = 0;
    size_t queryPosRemovedBack  = 0;
    size_t refPosRemovedFront   = 0;
    size_t refPosRemovedBack    = 0;

    // if mapped
    if (IsMapped()) {

        // update CIGAR - clip front ops, then clip back ops
        Cigar cigar = std::move(impl_.CigarData());
        size_t offsetRemaining = startOffset;
        while (offsetRemaining > 0 && !cigar.empty()) {
            CigarOperation& firstOp = cigar.front();
            const CigarOperationType firstOpType = firstOp.Type();
            const size_t firstOpLength = firstOp.Length();

            const bool shouldUpdateQueryPos = ((bam_cigar_type(static_cast<int>(firstOpType)) & 0x1) != 0);
            const bool shouldUpdateRefPos = ((bam_cigar_type(static_cast<int>(firstOpType)) & 0x2) != 0);

            if (firstOpLength <= offsetRemaining) {

                cigar.erase(cigar.begin());

                if (shouldUpdateQueryPos)
                    queryPosRemovedFront += firstOpLength;
                if (shouldUpdateRefPos)
                    refPosRemovedFront += firstOpLength;

                offsetRemaining -= firstOpLength;

            } else {

                firstOp.Length(firstOpLength - offsetRemaining);

                if (shouldUpdateQueryPos)
                    queryPosRemovedFront += offsetRemaining;
                if (shouldUpdateRefPos)
                    refPosRemovedFront += offsetRemaining;

                offsetRemaining = 0;
            }
        }

        offsetRemaining = endOffset;
        while (offsetRemaining > 0 && !cigar.empty()) {
            CigarOperation& lastOp = cigar.back();
            const CigarOperationType lastOpType = lastOp.Type();
            const size_t lastOpLength = lastOp.Length();

            const bool shouldUpdateQueryPos = ((bam_cigar_type(static_cast<int>(lastOpType)) & 0x1) != 0);
            const bool shouldUpdateRefPos = ((bam_cigar_type(static_cast<int>(lastOpType)) & 0x2) != 0);

            if (lastOpLength <= offsetRemaining) {
                cigar.pop_back();

                if (shouldUpdateQueryPos)
                    queryPosRemovedBack += lastOpLength;
                if (shouldUpdateRefPos)
                    refPosRemovedBack += lastOpLength;

                offsetRemaining -= lastOpLength;

            } else {
                lastOp.Length(lastOpLength - offsetRemaining);

                if (shouldUpdateQueryPos)
                    queryPosRemovedBack += offsetRemaining;
                if (shouldUpdateRefPos)
                    refPosRemovedBack += offsetRemaining;

                offsetRemaining = 0;
            }
        }
        impl_.CigarData(cigar);

        // update aligned reference position
        if (clipToQuery) {
            const Position origPosition = impl_.Position();
            impl_.Position(origPosition + refPosRemovedFront);
        } else {
            impl_.Position(start);
        }
    }

    const string origSequence = std::move(Sequence(Orientation::GENOMIC));
    const QualityValues origQualities = std::move(Qualities(Orientation::GENOMIC));

    size_t clipIndex;
    size_t clipLength;
    if (clipToQuery) {
        clipIndex = startOffset;
        clipLength = (end - start);
    } else {
        const size_t origSeqLength = origSequence.length();
        const size_t newSeqLength = (origSeqLength - queryPosRemovedBack) - queryPosRemovedFront;
        clipIndex = queryPosRemovedFront;
        clipLength = newSeqLength;
    }

    // clip seq, quals
    const string sequence = std::move(internal::Clip(origSequence, clipIndex, clipLength));
    const QualityValues qualities = std::move(internal::Clip(origQualities, clipIndex, clipLength));
    impl_.SetSequenceAndQualities(sequence, qualities.Fastq());

    // clip PacBio tags
    QualityValues deletionQV = std::move(internal::Clip(DeletionQV(Orientation::GENOMIC), clipIndex, clipLength));
    QualityValues insertionQV = std::move(internal::Clip(InsertionQV(Orientation::GENOMIC), clipIndex, clipLength));
    QualityValues mergeQV = std::move(internal::Clip(MergeQV(Orientation::GENOMIC), clipIndex, clipLength));
    QualityValues substitutionQV = std::move(internal::Clip(SubstitutionQV(Orientation::GENOMIC), clipIndex, clipLength));
    Frames ipd = std::move(internal::Clip(IPD(Orientation::GENOMIC).Data(), clipIndex, clipLength));
    Frames pulseWidth = std::move(internal::Clip(PulseWidth(Orientation::GENOMIC).Data(), clipIndex, clipLength));
    string deletionTag = std::move(internal::Clip(DeletionTag(Orientation::GENOMIC), clipIndex, clipLength));
    string substitutionTag = std::move(internal::Clip(SubstitutionTag(Orientation::GENOMIC), clipIndex, clipLength));

    // restore native orientation
    if (!isForwardStrand) {
        internal::Reverse(deletionQV);
        internal::Reverse(insertionQV);
        internal::Reverse(mergeQV);
        internal::Reverse(substitutionQV);
        internal::Reverse(ipd);
        internal::Reverse(pulseWidth);
        internal::ReverseComplement(deletionTag);
        internal::ReverseComplement(substitutionTag);
    }

    // update BAM tags
    TagCollection tags = impl_.Tags();
    tags[internal::tagName_deletionQV]      = deletionQV.Fastq();
    tags[internal::tagName_insertionQV]     = insertionQV.Fastq();
    tags[internal::tagName_mergeQV]         = mergeQV.Fastq();
    tags[internal::tagName_substitutionQV]  = substitutionQV.Fastq();
    tags[internal::tagName_ipd]             = ipd.Data();
    tags[internal::tagName_pulseWidth]      = pulseWidth.Data();
    tags[internal::tagName_deletionTag]     = deletionTag;
    tags[internal::tagName_substitutionTag] = substitutionTag;
    impl_.Tags(tags);

    // update query start/end
    if (clipToQuery) {
        internal::CreateOrEdit(internal::tagName_queryStart, start, &impl_);
        internal::CreateOrEdit(internal::tagName_queryEnd,   end,   &impl_);
    } else {
        if (isForwardStrand) {
            const Position qStart = origQStart + queryPosRemovedFront;
            const Position qEnd   = origQEnd   - queryPosRemovedBack;
            internal::CreateOrEdit(internal::tagName_queryStart, qStart, &impl_);
            internal::CreateOrEdit(internal::tagName_queryEnd,   qEnd,   &impl_);
        } else {
            const Position qStart = origQStart + queryPosRemovedBack;
            const Position qEnd   = origQEnd   - queryPosRemovedFront;
            internal::CreateOrEdit(internal::tagName_queryStart, qStart, &impl_);
            internal::CreateOrEdit(internal::tagName_queryEnd,   qEnd,   &impl_);
        }
    }

    // reset any cached aligned start/end
    alignedStart_ = PacBio::BAM::UnmappedPosition;
    alignedEnd_   = PacBio::BAM::UnmappedPosition;

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

string BamRecord::FetchBases(const string& tagName,
                             const Orientation orientation,
                             const bool aligned,
                             const bool exciseSoftClips) const
{
    const bool isBamSeq = (tagName == internal::tagName_SEQ);

    // fetch SAM/BAM SEQ field
    if (isBamSeq) {
        string seq = std::move(impl_.Sequence());

        // clip / gapify
        internal::MaybeClipAndGapifyBases(impl_,
                                          aligned,
                                          exciseSoftClips,
                                          seq);
        // rev-comp
        internal::MaybeReverseComplementSeq(isBamSeq,
                                            impl_.IsReverseStrand(),
                                            orientation,
                                            seq);
        return seq;
    }

    // other tags of 'bases' type
    else {

        const Tag& seqTag = impl_.TagValue(tagName);
        if (seqTag.IsNull())
            return string();

        bool ok;
        string seq = std::move(seqTag.ToString(&ok));
        if (!ok)
            return string();

        // rev-comp
        internal::MaybeReverseComplementSeq(isBamSeq,
                                            impl_.IsReverseStrand(),
                                            orientation,
                                            seq);
        // clip / gapify
        internal::MaybeClipAndGapifyBases(impl_,
                                          aligned,
                                          exciseSoftClips,
                                          seq);
        return seq;
    }
}

Frames BamRecord::FetchFrames(const string& tagName,
                              const Orientation orientation,
                              const bool aligned,
                              const bool exciseSoftClips) const
{
    Frames frames;
    const Tag& frameTag = impl_.TagValue(tagName);
    if (frameTag.IsNull())
        return frames;

    // lossy frame codes
    if (frameTag.IsUInt8Array()) {

        bool ok;
        const vector<uint8_t> codes = std::move(frameTag.ToUInt8Array(&ok));
        if (!ok)
            return Frames();
        frames = std::move(Frames::CodeToFrames(codes));
    }

    // lossless frame data
    else {
        assert(frameTag.IsUInt16Array());
        bool ok;
        const vector<uint16_t> losslessFrames = std::move(frameTag.ToUInt16Array(&ok));
        if (!ok)
            return Frames();
        frames.Data(std::move(losslessFrames));
    }

    // reverse, if needed
    internal::MaybeReverseFrames(impl_.IsReverseStrand(),
                                 orientation,
                                 frames);

    // clip / gapify
    internal::MaybeClipAndGapifyFrames(impl_,
                                      aligned,
                                      exciseSoftClips,
                                      frames);

    return frames;
}

QualityValues BamRecord::FetchQualities(const string& tagName,
                                        const Orientation orientation,
                                        const bool aligned,
                                        const bool exciseSoftClips) const
{
    const bool isBamQual = (tagName == internal::tagName_QUAL);

    // fetch SAM/BAM QUAL field
    if (isBamQual) {
        QualityValues quals = std::move(impl_.Qualities());

        // clip / gapify
        internal::MaybeClipAndGapifyQualities(impl_,
                                              aligned,
                                              exciseSoftClips,
                                              quals);
        // rev-comp
        internal::MaybeReverseQuals(isBamQual,
                                    impl_.IsReverseStrand(),
                                    orientation,
                                    quals);
        return quals;
    }

    // other tags of 'qualities' type
    else {

        // fetch data
        const Tag& qvsTag = impl_.TagValue(tagName);
        if (qvsTag.IsNull())
            return QualityValues();
        bool ok;
        QualityValues quals = std::move(QualityValues::FromFastq(qvsTag.ToString(&ok)));
        if (!ok)
            return QualityValues();

        // rev-comp
        internal::MaybeReverseQuals(isBamQual,
                                    impl_.IsReverseStrand(),
                                    orientation,
                                    quals);
        // clip / gapify
        internal::MaybeClipAndGapifyQualities(impl_,
                                              aligned,
                                              exciseSoftClips,
                                              quals);
        return quals;
    }
}

string BamRecord::FullName(void) const
{ return impl_.Name(); }

bool BamRecord::HasDeletionQV(void) const
{ return impl_.HasTag(internal::tagName_deletionQV); }

bool BamRecord::HasDeletionTag(void) const
{ return impl_.HasTag(internal::tagName_deletionTag); }

bool BamRecord::HasInsertionQV(void) const
{ return impl_.HasTag(internal::tagName_insertionQV); }

bool BamRecord::HasIPD(void) const
{ return impl_.HasTag(internal::tagName_ipd); }

bool BamRecord::HasMergeQV(void) const
{ return impl_.HasTag(internal::tagName_mergeQV); }

bool BamRecord::HasPulseWidth(void) const
{ return impl_.HasTag(internal::tagName_pulseWidth); }

bool BamRecord::HasSubstitutionQV(void) const
{ return impl_.HasTag(internal::tagName_substitutionQV); }

bool BamRecord::HasSubstitutionTag(void) const
{ return impl_.HasTag(internal::tagName_substitutionTag); }

BamHeader::SharedPtr BamRecord::Header(void) const
{ return header_; }

int32_t BamRecord::HoleNumber(void) const
{
    const Tag& holeNumber = impl_.TagValue(internal::tagName_holeNumber);
    bool ok;
    int32_t result = holeNumber.ToInt32(&ok);
    if (!ok)
        return -1;
    return result;
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
        internal::CreateOrEdit(internal::tagName_ipd, frames.Downsampled(), &impl_);
    else
        internal::CreateOrEdit(internal::tagName_ipd, frames.Data(), &impl_);
    return *this;
}

bool BamRecord::IsMapped(void) const
{ return impl_.IsMapped(); }

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

string BamRecord::MovieName(void) const
{ return ReadGroup().MovieName(); }

int32_t BamRecord::NumPasses(void) const
{
    const Tag& numPasses = impl_.TagValue(internal::tagName_numPasses);
    bool ok;
    const int32_t result = numPasses.ToInt32(&ok);
    if (!ok)
        return -1;
    return result;
}

BamRecord& BamRecord::NumPasses(const int32_t numPasses)
{
    internal::CreateOrEdit(internal::tagName_numPasses, numPasses, &impl_);
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
        internal::CreateOrEdit(internal::tagName_pulseWidth, frames.Downsampled(), &impl_);
    else
        internal::CreateOrEdit(internal::tagName_pulseWidth, frames.Data(), &impl_);
    return *this;
}

Position BamRecord::QueryEnd(void) const
{
    const Tag& qe = impl_.TagValue(internal::tagName_queryEnd);
    if (!qe.IsNull()) {
        bool ok;
        Position result = qe.ToInt32(&ok);
        if (!ok)
            return PacBio::BAM::UnmappedPosition;
        return result;
    }

    // missing qe tag - check movie name

    // otherwise fallback to SAM/BAM::POS + length
    return impl_.Position() + impl_.Sequence().length();
}

Position BamRecord::QueryStart(void) const
{
    const Tag& qs = impl_.TagValue(internal::tagName_queryStart);
    if (!qs.IsNull()) {
        bool ok;
        Position result = qs.ToInt32(&ok);
        if (!ok)
            return PacBio::BAM::UnmappedPosition;
        return result;
    }

    // missing qs tag - chcek movie name

    // otherwise fallback to SAM/BAM::POS
    return impl_.Position();
}

Accuracy BamRecord::ReadAccuracy(void) const
{
    const Tag& readAccuracy = impl_.TagValue(internal::tagName_readAccuracy);
    return Accuracy(readAccuracy.ToInt32());
}

BamRecord& BamRecord::ReadAccuracy(const Accuracy& accuracy)
{
    internal::CreateOrEdit(internal::tagName_readAccuracy,
                           (int32_t)accuracy,
                           &impl_);
    return *this;
}

ReadGroupInfo BamRecord::ReadGroup(void) const
{
    if (!header_)
        return ReadGroupInfo();
    return header_->ReadGroup(ReadGroupId());
}

string BamRecord::ReadGroupId(void) const
{
    const Tag& rgTag = impl_.TagValue(internal::tagName_readGroup);
    return rgTag.ToString();
}

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

Position BamRecord::ReferenceStart(void) const
{ return impl_.Position(); }

std::string BamRecord::Sequence(const Orientation orientation,
                                bool aligned,
                                bool exciseSoftClips) const
{
    return FetchBases("SEQ",
                      orientation,
                      aligned,
                      exciseSoftClips);
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
    const string& typeName = ReadGroup().ReadType();
    return internal::NameToType(typeName);
}



//BamRecord& BamRecord::MovieName(const string& movieName)
//{
//    UpdateName();
//    return *this;
//}

//BamRecord& BamRecord::QueryEnd(const Position pos)
//{
//    internal::CreateOrEdit(internal::tagName_queryEnd,
//                           (int32_t)pos,
//                           &impl_);
//    UpdateName();
//    return *this;
//}

//BamRecord& BamRecord::QueryStart(const Position pos)
//{
//    internal::CreateOrEdit(internal::tagName_queryStart,
//                           (int32_t)pos,
//                           &impl_);
//    UpdateName();
//    return *this;
//}

//BamRecord& BamRecord::ReadGroup(const ReadGroupInfo& rg)
//{
//    internal::CreateOrEdit(internal::tagName_readGroup, rg.Id(), &impl_);
//    // update anything else?
//    return *this;
//}

//BamRecord& BamRecord::ReadGroupId(const std::string& id)
//{
//    internal::CreateOrEdit(internal::tagName_readGroup,
//                           id,
//                           &impl_);
//    return *this;
//}

//BamRecord& BamRecord::ReferenceStart(const Position pos)
//{
//    impl_.Position(pos);
//    return *this;
//}

