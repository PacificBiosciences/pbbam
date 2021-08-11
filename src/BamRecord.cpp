#include "PbbamInternalConfig.h"

#include <pbbam/BamRecord.h>

#include <cassert>
#include <cstddef>
#include <cstdint>

#include <stdexcept>

#include <htslib/sam.h>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include <pbcopper/data/Clipping.h>
#include <pbcopper/data/FrameEncoders.h>
#include <pbcopper/data/internal/ClippingImpl.h>

#include <pbbam/StringUtilities.h>
#include <pbbam/ZmwTypeMap.h>
#include <pbbam/virtual/VirtualRegionTypeMap.h>

#include "BamRecordTags.h"
#include "MemoryUtils.h"
#include "Pulse2BaseCache.h"
#include "SequenceUtils.h"

namespace PacBio {
namespace BAM {
namespace {

// record type names
const std::string recordTypeName_ZMW{"ZMW"};
const std::string recordTypeName_Polymerase{"POLYMERASE"};
const std::string recordTypeName_HqRegion{"HQREGION"};
const std::string recordTypeName_Subread{"SUBREAD"};
const std::string recordTypeName_CCS{"CCS"};
const std::string recordTypeName_Scrap{"SCRAP"};
const std::string recordTypeName_Transcript{"TRANSCRIPT"};
const std::string recordTypeName_Unknown{"UNKNOWN"};

int32_t HoleNumberFromName(const std::string& fullName)
{
    const auto mainTokens = Split(fullName, '/');
    if (mainTokens.size() < 2) {
        throw std::runtime_error{"[pbbam] BAM record ERROR: malformed record name: " + fullName};
    }

    try {
        return std::stoi(mainTokens.at(1));
    } catch (const std::exception&) {
        throw std::runtime_error{"[pbbam] BAM record ERROR: invalid hole number: '" +
                                 mainTokens.at(1) + "'"};
    }
}

Data::Position QueryEndFromName(const std::string& fullName)
{
    const auto mainTokens = Split(fullName, '/');
    if (mainTokens.size() != 3) {
        throw std::runtime_error{"[pbbam] BAM record ERROR: malformed record name: " + fullName};
    }

    const auto queryTokens = Split(mainTokens.at(2), '_');
    if (queryTokens.size() != 2) {
        throw std::runtime_error{"[pbbam] BAM record ERROR: malformed record name: " + fullName};
    }

    return std::stoi(queryTokens.at(1));
}

Data::Position QueryStartFromName(const std::string& fullName)
{
    const auto mainTokens = Split(fullName, '/');
    if (mainTokens.size() != 3) {
        throw std::runtime_error{"[pbbam] BAM record ERROR: malformed record name: " + fullName};
    }

    const auto queryTokens = Split(mainTokens.at(2), '_');
    if (queryTokens.size() != 2) {
        throw std::runtime_error{"[pbbam] BAM record ERROR: malformed record name: " + fullName};
    }

    return std::stoi(queryTokens.at(0));
}

std::string Label(const BamRecordTag tag) { return BamRecordTags::LabelFor(tag); }

BamRecordImpl* CreateOrEdit(const BamRecordTag tag, const Tag& value, BamRecordImpl* impl)
{
    if (impl->HasTag(tag)) {
        impl->EditTag(tag, value);
    } else {
        impl->AddTag(tag, value);
    }
    return impl;
}

std::pair<int32_t, int32_t> AlignedOffsets(const BamRecord& record, const int seqLength)
{
    int32_t startOffset = 0;
    int32_t endOffset = seqLength;

    const auto& b = BamRecordMemory::GetRawData(record);
    uint32_t* cigarData = bam_get_cigar(b.get());
    const size_t numCigarOps = b->core.n_cigar;
    if (numCigarOps > 0) {

        // start offset
        for (size_t i = 0; i < numCigarOps; ++i) {
            const auto type = static_cast<Data::CigarOperationType>(bam_cigar_op(cigarData[i]));
            if (type == Data::CigarOperationType::HARD_CLIP) {
                if (startOffset != 0 && startOffset != seqLength) {
                    startOffset = -1;
                    break;
                }
            } else if (type == Data::CigarOperationType::SOFT_CLIP) {
                startOffset += bam_cigar_oplen(cigarData[i]);
            } else {
                break;
            }
        }

        // end offset
        for (int i = numCigarOps - 1; i >= 0; --i) {
            const auto type = static_cast<Data::CigarOperationType>(bam_cigar_op(cigarData[i]));
            if (type == Data::CigarOperationType::HARD_CLIP) {
                if (endOffset != 0 && endOffset != seqLength) {
                    endOffset = -1;
                    break;
                }
            } else if (type == Data::CigarOperationType::SOFT_CLIP) {
                endOffset -= bam_cigar_oplen(cigarData[i]);
            } else {
                break;
            }
        }

        if (endOffset == 0) {
            endOffset = seqLength;
        }
    }
    return {startOffset, endOffset};
}

template <class InputIt, class Size, class OutputIt>
OutputIt Move_N(InputIt first, Size count, OutputIt result)
{
    return std::move(first, first + count, result);
}

template <typename T>
T ClipSeqQV(const T& input, const size_t pos, const size_t len)
{
    if (input.empty()) {
        return {};
    }
    return T{input.cbegin() + pos, input.cbegin() + pos + len};
}

template <typename T>
T ClipPulse(const T& input, Pulse2BaseCache* p2bCache, const size_t pos, const size_t len)
{
    assert(p2bCache);
    if (input.empty()) {
        return {};
    }

    // find start
    size_t start = p2bCache->FindFirst();
    size_t basesSeen = 0;
    while (basesSeen < pos) {
        start = p2bCache->FindNext(start);
        ++basesSeen;
    }

    // find end
    size_t end = start;
    size_t seen = 1;
    while (seen < len) {
        end = p2bCache->FindNext(end);
        ++seen;
    }

    // return clipped
    return {input.cbegin() + start, input.cbegin() + end + 1};
}

template <typename F, typename N>
void ClipAndGapify(const BamRecordImpl& impl, const bool aligned, const bool exciseSoftClips,
                   F* seq, N paddingNullValue, N deletionNullValue)
{
    assert(seq);

    const bool clipOrGapRequested = aligned || exciseSoftClips;
    if (impl.IsMapped() && clipOrGapRequested) {
        // determine final container length
        auto incrementsOutputLength = [](const Data::CigarOperationType type, const bool isAligned,
                                         const bool exciseSoftClipsFromAln) {
            if (type == Data::CigarOperationType::HARD_CLIP ||
                type == Data::CigarOperationType::REFERENCE_SKIP) {
                return false;
            } else if (type == Data::CigarOperationType::SOFT_CLIP && exciseSoftClipsFromAln) {
                return false;
            } else if (!isAligned && (type == Data::CigarOperationType::DELETION ||
                                      type == Data::CigarOperationType::PADDING)) {
                return false;
            } else {
                return true;
            }
        };

        size_t outputLength = 0;
        const auto cigar = impl.CigarData();
        for (const auto& op : cigar) {
            if (incrementsOutputLength(op.Type(), aligned, exciseSoftClips)) {
                outputLength += op.Length();
            }
        }

        // move original data to temp, prep output container size
        F originalSeq = std::move(*seq);
        seq->resize(outputLength);

        // apply CIGAR ops
        size_t srcIndex = 0;
        size_t dstIndex = 0;
        for (const auto& op : cigar) {
            const auto opType = op.Type();
            const auto opLength = op.Length();

            // nothing to do for hard-clipped & ref-skipped positions
            if (opType == Data::CigarOperationType::HARD_CLIP ||
                opType == Data::CigarOperationType::REFERENCE_SKIP) {
                continue;
            }

            // maybe skip soft-clipped positions
            else if (opType == Data::CigarOperationType::SOFT_CLIP) {
                if (exciseSoftClips) {
                    srcIndex += opLength;
                } else {
                    Move_N(originalSeq.begin() + srcIndex, opLength, seq->begin() + dstIndex);
                    srcIndex += opLength;
                    dstIndex += opLength;
                }
            }

            // maybe add deletion/padding values
            // either way, srcIndex is not incremented
            else if (opType == Data::CigarOperationType::DELETION) {
                if (aligned) {
                    for (size_t i = 0; i < opLength; ++i) {
                        (*seq)[dstIndex] = deletionNullValue;
                        ++dstIndex;
                    }
                }
            } else if (opType == Data::CigarOperationType::PADDING) {
                if (aligned) {
                    for (size_t i = 0; i < opLength; ++i) {
                        (*seq)[dstIndex] = paddingNullValue;
                        ++dstIndex;
                    }
                }
            }

            // all other CIGAR ops
            else {
                Move_N(originalSeq.begin() + srcIndex, opLength, seq->begin() + dstIndex);
                srcIndex += opLength;
                dstIndex += opLength;
            }
        }
    }
}

void ClipAndGapifyBases(const BamRecordImpl& impl, const bool aligned, const bool exciseSoftClips,
                        std::string* seq)
{
    ClipAndGapify<std::string, char>(impl, aligned, exciseSoftClips, seq, '*', '-');
}

void ClipAndGapifyFrames(const BamRecordImpl& impl, const bool aligned, const bool exciseSoftClips,
                         Data::Frames* frames)
{
    assert(frames);
    std::vector<uint16_t> data{std::move(frames->Data())};
    ClipAndGapify<std::vector<uint16_t>, uint16_t>(impl, aligned, exciseSoftClips, &data, 0, 0);
    frames->Data(data);
}

void ClipAndGapifyPhotons(const BamRecordImpl& impl, const bool aligned, const bool exciseSoftClips,
                          std::vector<float>* data)
{
    ClipAndGapify<std::vector<float>, float>(impl, aligned, exciseSoftClips, data, 0.0, 0.0);
}

void ClipAndGapifyQualities(const BamRecordImpl& impl, const bool aligned,
                            const bool exciseSoftClips, Data::QualityValues* quals)
{
    ClipAndGapify<Data::QualityValues, Data::QualityValue>(
        impl, aligned, exciseSoftClips, quals, Data::QualityValue{0}, Data::QualityValue{0});
}

void ClipAndGapifyUInts(const BamRecordImpl& impl, const bool aligned, const bool exciseSoftClips,
                        std::vector<uint32_t>* data)
{
    ClipAndGapify<std::vector<uint32_t>, uint32_t>(impl, aligned, exciseSoftClips, data, 0, 0);
}

void ClipAndGapifyUInt8s(const BamRecordImpl& impl, const bool aligned, const bool exciseSoftClips,
                         std::vector<uint8_t>* data)
{
    ClipAndGapify<std::vector<uint8_t>, uint8_t>(impl, aligned, exciseSoftClips, data, 0, 0);
}

Data::FrameEncoder IpdEncoder(const BamRecord& record)
{
    try {
        return record.ReadGroup().IpdFrameEncoder();
    } catch (const std::exception&) {
        // fallback to V1 in corner cases w/ no read group set
        return Data::V1FrameEncoder{};
    }
}

Data::FrameEncoder PwEncoder(const BamRecord& record)
{
    try {
        return record.ReadGroup().PulseWidthFrameEncoder();
    } catch (const std::exception&) {
        // fallback to V1 in corner cases w/ no read group set
        return Data::V1FrameEncoder{};
    }
}

RecordType NameToType(const std::string& name)
{
    if (name == recordTypeName_Subread) {
        return RecordType::SUBREAD;
    }
    if (name == recordTypeName_ZMW || name == recordTypeName_Polymerase) {
        return RecordType::ZMW;
    }
    if (name == recordTypeName_HqRegion) {
        return RecordType::HQREGION;
    }
    if (name == recordTypeName_CCS) {
        return RecordType::CCS;
    }
    if (name == recordTypeName_Scrap) {
        return RecordType::SCRAP;
    }
    if (name == recordTypeName_Transcript) {
        return RecordType::TRANSCRIPT;
    }
    return RecordType::UNKNOWN;
}

void OrientBasesAsRequested(std::string* bases, Data::Orientation current,
                            Data::Orientation requested, bool isReverseStrand, bool isPulse)
{
    assert(bases);
    if (current != requested && isReverseStrand) {
        if (isPulse) {
            ReverseComplementCaseSens(*bases);
        } else {
            ReverseComplement(*bases);
        }
    }
}

template <typename Container>
void OrientTagDataAsRequested(Container* data, Data::Orientation current,
                              Data::Orientation requested, bool isReverseStrand)
{
    assert(data);
    if (current != requested && isReverseStrand) {
        std::reverse(data->begin(), data->end());
    }
}

}  // namespace

const float BamRecord::photonFactor = 10.0;

BamRecord::BamRecord() = default;

BamRecord::BamRecord(BamHeader header) : header_{std::move(header)} {}

BamRecord::BamRecord(BamRecordImpl impl) : impl_{std::move(impl)} {}

BamRecord::BamRecord(const BamRecord& other)
    : impl_{other.impl_}
    , header_{other.header_}
    , alignedStart_{other.alignedStart_}
    , alignedEnd_{other.alignedEnd_}
{
}

BamRecord::BamRecord(BamRecord&&) noexcept = default;

BamRecord& BamRecord::operator=(const BamRecord& other)
{
    if (this != &other) {
        impl_ = other.impl_;
        header_ = other.header_;
        alignedStart_ = other.alignedStart_;
        alignedEnd_ = other.alignedEnd_;
        p2bCache_.reset();  // just reset, for now at least
    }
    return *this;
}

BamRecord& BamRecord::operator=(BamRecord&&) noexcept = default;

BamRecord::~BamRecord() = default;

Data::Position BamRecord::AlignedEnd() const
{
    if (alignedEnd_ == Data::UnmappedPosition) {
        CalculateAlignedPositions();
    }
    return alignedEnd_;
}

Data::Position BamRecord::AlignedStart() const
{
    if (alignedStart_ == Data::UnmappedPosition) {
        CalculateAlignedPositions();
    }
    return alignedStart_;
}

Data::Strand BamRecord::AlignedStrand() const
{
    return impl_.IsReverseStrand() ? Data::Strand::REVERSE : Data::Strand::FORWARD;
}

Data::QualityValues BamRecord::AltLabelQV(Data::Orientation orientation, bool aligned,
                                          bool exciseSoftClips, PulseBehavior pulseBehavior) const
{
    return FetchQualities(BamRecordTag::ALT_LABEL_QV, orientation, aligned, exciseSoftClips,
                          pulseBehavior);
}

BamRecord& BamRecord::AltLabelQV(const Data::QualityValues& altLabelQVs)
{
    CreateOrEdit(BamRecordTag::ALT_LABEL_QV, altLabelQVs.Fastq(), &impl_);
    return *this;
}

std::string BamRecord::AltLabelTag(Data::Orientation orientation, bool aligned,
                                   bool exciseSoftClips, PulseBehavior pulseBehavior) const
{
    return FetchBases(BamRecordTag::ALT_LABEL_TAG, orientation, aligned, exciseSoftClips,
                      pulseBehavior);
}

BamRecord& BamRecord::AltLabelTag(const std::string& tags)
{
    CreateOrEdit(BamRecordTag::ALT_LABEL_TAG, tags, &impl_);
    return *this;
}

int16_t BamRecord::BarcodeForward() const { return Barcodes().first; }

int16_t BamRecord::BarcodeReverse() const { return Barcodes().second; }

uint8_t BamRecord::BarcodeQuality() const
{
    const auto tagName = BamRecordTags::LabelFor(BamRecordTag::BARCODE_QUALITY);
    const auto bq = impl_.TagValue(tagName);
    if (bq.IsNull()) {
        return 0;  // ?? "missing" value for tags ?? should we consider boost::optional<T> for these kind of guys ??
    }
    return bq.ToUInt8();
}

BamRecord& BamRecord::BarcodeQuality(const uint8_t quality)
{
    CreateOrEdit(BamRecordTag::BARCODE_QUALITY, quality, &impl_);
    return *this;
}

std::pair<int16_t, int16_t> BamRecord::Barcodes() const
{
    const auto tagName = BamRecordTags::LabelFor(BamRecordTag::BARCODES);
    const auto bc = impl_.TagValue(tagName);
    if (bc.IsNull()) {
        throw std::runtime_error{
            "[pbbam] BAM record ERROR: barcode tag (bc) was requested but is missing"};
    }

    // NOTE: barcodes are still stored, per the spec, as uint16, even though
    // we're now using them as int16_t in the API (bug 31511)
    //
    if (!bc.IsUInt16Array()) {
        throw std::runtime_error{
            "[pbbam] BAM record ERROR: barcode tag (bc) is malformed: should be a uint16_t array "
            "of size==2."};
    }
    const auto bcArray = bc.ToUInt16Array();
    if (bcArray.size() != 2) {
        throw std::runtime_error{
            "[pbbam] BAM record ERROR: barcode tag (bc) is malformed: should be a uint16_t array "
            "of size==2."};
    }

    return {boost::numeric_cast<int16_t>(bcArray[0]), boost::numeric_cast<int16_t>(bcArray[1])};
}

BamRecord& BamRecord::Barcodes(const std::pair<int16_t, int16_t>& barcodeIds)
{
    const std::vector<uint16_t> data{boost::numeric_cast<uint16_t>(barcodeIds.first),
                                     boost::numeric_cast<uint16_t>(barcodeIds.second)};
    CreateOrEdit(BamRecordTag::BARCODES, data, &impl_);
    return *this;
}

void BamRecord::CalculateAlignedPositions() const
{
    // reset
    ResetCachedPositions();

    // skip if unmapped, or has no queryStart/End
    if (!impl_.IsMapped()) {
        return;
    }

    // get the query start/end
    const auto seqLength = static_cast<int>(impl_.SequenceLength());
    const bool isCcsOrTranscript = IsCcsOrTranscript(Type());
    const Data::Position qStart = isCcsOrTranscript ? 0 : QueryStart();
    const Data::Position qEnd = isCcsOrTranscript ? seqLength : QueryEnd();

    if (qStart == Data::UnmappedPosition || qEnd == Data::UnmappedPosition) {
        return;
    }

    // determine clipped end ranges
    const auto alignedOffsets = AlignedOffsets(*this, seqLength);
    const auto startOffset = alignedOffsets.first;
    const auto endOffset = alignedOffsets.second;
    if (endOffset == -1 || startOffset == -1) {
        return;  // TODO: handle error more??
    }

    // store aligned positions (polymerase read coordinates)
    if (impl_.IsReverseStrand()) {
        alignedStart_ = qStart + (seqLength - endOffset);
        alignedEnd_ = qEnd - startOffset;
    } else {
        alignedStart_ = qStart + startOffset;
        alignedEnd_ = qEnd - (seqLength - endOffset);
    }
}

void BamRecord::CalculatePulse2BaseCache() const
{
    // skip already calculated
    if (p2bCache_) {
        return;
    }

    // else try to calculate p2b cache.
    if (!HasPulseCall()) {
        throw std::runtime_error{
            "[pbbam] BAM record ERROR: cannot calculate pulse2base mapping without 'pc' tag."};
    }
    const auto pulseCalls = FetchBases(BamRecordTag::PULSE_CALL, Data::Orientation::NATIVE, false,
                                       false, PulseBehavior::ALL);
    p2bCache_ = std::make_unique<Pulse2BaseCache>(pulseCalls);
}

Data::Cigar BamRecord::CigarData(bool exciseAllClips) const
{
    auto isClippingOp = [](const auto& op) {
        const auto type = op.Type();
        return type == Data::CigarOperationType::SOFT_CLIP ||
               type == Data::CigarOperationType::HARD_CLIP;
    };

    auto cigar = impl_.CigarData();
    if (exciseAllClips) {
        cigar.erase(std::remove_if(cigar.begin(), cigar.end(), isClippingOp), cigar.end());
    }
    return cigar;
}

BamRecord& BamRecord::Clip(const ClipType clipType, const Data::Position start,
                           const Data::Position end, const bool exciseFlankingInserts)
{
    switch (clipType) {
        case ClipType::CLIP_NONE:
            return *this;
        case ClipType::CLIP_TO_QUERY:
            // exciseFlankingInserts ignored, just clipping to query coordinates
            return ClipToQuery(start, end);
        case ClipType::CLIP_TO_REFERENCE:
            return ClipToReference(start, end, exciseFlankingInserts);
        default:
            throw std::runtime_error{"[pbbam] BAM record ERROR: unsupported clip type requested"};
    }
}

BamRecord BamRecord::Clipped(const BamRecord& input, const ClipType clipType,
                             const Data::Position start, const Data::Position end,
                             const bool exciseFlankingInserts)
{
    return input.Clipped(clipType, start, end, exciseFlankingInserts);
}

BamRecord BamRecord::Clipped(const ClipType clipType, const Data::Position start,
                             const Data::Position end, const bool exciseFlankingInserts) const
{
    BamRecord result(*this);
    result.Clip(clipType, start, end, exciseFlankingInserts);
    return result;
}

void BamRecord::ClipTags(const size_t clipFrom, const size_t clipLength)
{
    const auto rg = ReadGroup();
    const auto ipdCodec = rg.IpdCodec();
    const auto ipdEncoder = rg.IpdFrameEncoder();
    const auto pwCodec = rg.PulseWidthCodec();
    const auto pwEncoder = rg.IpdFrameEncoder();

    // update BAM tags
    TagCollection tags = impl_.Tags();
    if (HasDeletionQV()) {
        tags[Label(BamRecordTag::DELETION_QV)] =
            ClipSeqQV(DeletionQV(Data::Orientation::NATIVE), clipFrom, clipLength).Fastq();
    }
    if (HasInsertionQV()) {
        tags[Label(BamRecordTag::INSERTION_QV)] =
            ClipSeqQV(InsertionQV(Data::Orientation::NATIVE), clipFrom, clipLength).Fastq();
    }
    if (HasMergeQV()) {
        tags[Label(BamRecordTag::MERGE_QV)] =
            ClipSeqQV(MergeQV(Data::Orientation::NATIVE), clipFrom, clipLength).Fastq();
    }
    if (HasSubstitutionQV()) {
        tags[Label(BamRecordTag::SUBSTITUTION_QV)] =
            ClipSeqQV(SubstitutionQV(Data::Orientation::NATIVE), clipFrom, clipLength).Fastq();
    }
    if (HasIPD()) {
        const auto label = Label(BamRecordTag::IPD);
        const auto ipd = IPD(Data::Orientation::NATIVE).Data();
        if (ipdCodec == Data::FrameCodec::RAW) {
            tags[label] = ClipSeqQV(ipd, clipFrom, clipLength);
        } else {
            tags[label] = ClipSeqQV(ipdEncoder.Encode(ipd), clipFrom, clipLength);
        }
    }
    if (HasPulseWidth()) {
        const auto label = Label(BamRecordTag::PULSE_WIDTH);
        const auto pw = PulseWidth(Data::Orientation::NATIVE).Data();
        if (pwCodec == Data::FrameCodec::RAW) {
            tags[label] = ClipSeqQV(pw, clipFrom, clipLength);
        } else {
            tags[label] = ClipSeqQV(pwEncoder.Encode(pw), clipFrom, clipLength);
        }
    }
    if (HasDeletionTag()) {
        tags[Label(BamRecordTag::DELETION_TAG)] =
            ClipSeqQV(DeletionTag(Data::Orientation::NATIVE), clipFrom, clipLength);
    }
    if (HasSubstitutionTag()) {
        tags[Label(BamRecordTag::SUBSTITUTION_TAG)] =
            ClipSeqQV(SubstitutionTag(Data::Orientation::NATIVE), clipFrom, clipLength);
    }

    // internal BAM tags
    if (HasPulseCall()) {

        // ensure p2bCache initialized
        CalculatePulse2BaseCache();
        Pulse2BaseCache* p2bCache = p2bCache_.get();

        if (HasAltLabelQV()) {
            tags[Label(BamRecordTag::ALT_LABEL_QV)] =
                ClipPulse(AltLabelQV(Data::Orientation::NATIVE), p2bCache, clipFrom, clipLength)
                    .Fastq();
        }
        if (HasLabelQV()) {
            tags[Label(BamRecordTag::LABEL_QV)] =
                ClipPulse(LabelQV(Data::Orientation::NATIVE), p2bCache, clipFrom, clipLength)
                    .Fastq();
        }
        if (HasPulseMergeQV()) {
            tags[Label(BamRecordTag::PULSE_MERGE_QV)] =
                ClipPulse(PulseMergeQV(Data::Orientation::NATIVE), p2bCache, clipFrom, clipLength)
                    .Fastq();
        }
        if (HasAltLabelTag()) {
            tags[Label(BamRecordTag::ALT_LABEL_TAG)] =
                ClipPulse(AltLabelTag(Data::Orientation::NATIVE), p2bCache, clipFrom, clipLength);
        }
        if (HasPulseCall()) {
            tags[Label(BamRecordTag::PULSE_CALL)] =
                ClipPulse(PulseCall(Data::Orientation::NATIVE), p2bCache, clipFrom, clipLength);
        }
        if (HasPkmean()) {
            tags[Label(BamRecordTag::PKMEAN)] = EncodePhotons(
                ClipPulse(Pkmean(Data::Orientation::NATIVE), p2bCache, clipFrom, clipLength));
        }
        if (HasPkmid()) {
            tags[Label(BamRecordTag::PKMID)] = EncodePhotons(
                ClipPulse(Pkmid(Data::Orientation::NATIVE), p2bCache, clipFrom, clipLength));
        }
        if (HasPkmean2()) {
            tags[Label(BamRecordTag::PKMEAN_2)] = EncodePhotons(
                ClipPulse(Pkmean2(Data::Orientation::NATIVE), p2bCache, clipFrom, clipLength));
        }
        if (HasPkmid2()) {
            tags[Label(BamRecordTag::PKMID_2)] = EncodePhotons(
                ClipPulse(Pkmid2(Data::Orientation::NATIVE), p2bCache, clipFrom, clipLength));
        }
        if (HasPrePulseFrames()) {
            tags[Label(BamRecordTag::PRE_PULSE_FRAMES)] = ClipPulse(
                PrePulseFrames(Data::Orientation::NATIVE).Data(), p2bCache, clipFrom, clipLength);
        }
        if (HasPulseCallWidth()) {
            tags[Label(BamRecordTag::PULSE_CALL_WIDTH)] = ClipPulse(
                PulseCallWidth(Data::Orientation::NATIVE).Data(), p2bCache, clipFrom, clipLength);
        }
        if (HasStartFrame()) {
            tags[Label(BamRecordTag::START_FRAME)] =
                ClipPulse(StartFrame(Data::Orientation::NATIVE), p2bCache, clipFrom, clipLength);
        }
    }

    impl_.Tags(tags);
}

void BamRecord::ClipFields(const size_t clipFrom, const size_t clipLength)
{
    const bool isForwardStrand = (AlignedStrand() == Data::Strand::FORWARD);

    // clip seq, quals
    std::string sequence{ClipSeqQV(Sequence(Data::Orientation::NATIVE), clipFrom, clipLength)};
    Data::QualityValues qualities{
        ClipSeqQV(Qualities(Data::Orientation::NATIVE), clipFrom, clipLength)};
    if (!isForwardStrand) {
        ReverseComplement(sequence);
        Reverse(qualities);
    }
    impl_.SetSequenceAndQualities(sequence, qualities.Fastq());

    ClipTags(clipFrom, clipLength);
}

BamRecord& BamRecord::ClipToQuery(const Data::Position start, const Data::Position end)
{
    // cache original coords, skip out if clip not needed
    const size_t seqLength = impl_.SequenceLength();
    const bool isCcsOrTranscript = IsCcsOrTranscript(Type());
    const Data::Position origQStart = isCcsOrTranscript ? 0 : QueryStart();
    const Data::Position origQEnd = isCcsOrTranscript ? seqLength : QueryEnd();
    if (start <= origQStart && end >= origQEnd) {
        return *this;
    }

    // calculate clipping
    Data::ClipToQueryConfig clipConfig{
        impl_.SequenceLength(), origQStart,      origQEnd,          start,           end,
        impl_.Position(),       AlignedStrand(), impl_.CigarData(), impl_.IsMapped()};
    auto result = Data::ClipToQuery(clipConfig);

    // update alignment info
    if (IsMapped()) {
        impl_.CigarData(std::move(result.cigar_));
        impl_.Position(result.refPos_);
    }

    // clip SEQ, QUAL, tags
    const auto clipFrom = result.clipOffset_;
    const auto clipLength = (end - start);
    ClipFields(clipFrom, clipLength);

    // update query start/end
    // TODO: update name to reflect new QS/QE ???
    CreateOrEdit(BamRecordTag::QUERY_START, start, &impl_);
    CreateOrEdit(BamRecordTag::QUERY_END, end, &impl_);

    // reset any cached aligned start/end
    ResetCachedPositions();
    return *this;
}

BamRecord& BamRecord::ClipToReference(const Data::Position start, const Data::Position end,
                                      const bool exciseFlankingInserts)
{
    // skip if not mapped, clipping to reference doesn't make sense
    // or should we even consider throwing here?
    if (!IsMapped()) {
        return *this;
    }

    // cache original coords
    const int seqLength = static_cast<int>(impl_.SequenceLength());
    const bool isCcsOrTranscript = IsCcsOrTranscript(Type());
    const Data::Position origQStart = isCcsOrTranscript ? 0 : QueryStart();
    const Data::Position origQEnd = isCcsOrTranscript ? seqLength : QueryEnd();
    const Data::Position origTStart = ReferenceStart();
    const Data::Position origTEnd = ReferenceEnd();

    // skip if already within requested clip range
    if (start <= origTStart && end >= origTEnd) {
        return *this;
    }
    assert(AlignedStart() >= origQStart);
    assert(AlignedEnd() <= origQEnd);

    // calculate clipping
    Data::ClipToReferenceConfig clipConfig{
        Data::ClipToQueryConfig{impl_.SequenceLength(), origQStart, origQEnd, start, end,
                                impl_.Position(), AlignedStrand(), impl_.CigarData(),
                                impl_.IsMapped()},
        ReferenceEnd(), start, end, exciseFlankingInserts};
    auto result = Data::ClipToReference(clipConfig);

    // update CIGAR and position
    impl_.CigarData(std::move(result.cigar_));
    impl_.Position(result.refPos_);

    // clip SEQ, QUAL, tags
    const Data::Position qStart = result.qStart_;
    const Data::Position qEnd = result.qEnd_;
    const size_t clipFrom = result.clipOffset_;
    const size_t clipLength = qEnd - qStart;
    ClipFields(clipFrom, clipLength);

    // update query start/end
    CreateOrEdit(BamRecordTag::QUERY_START, qStart, &impl_);
    CreateOrEdit(BamRecordTag::QUERY_END, qEnd, &impl_);

    // reset any cached aligned start/end
    ResetCachedPositions();
    return *this;
}

Data::QualityValues BamRecord::DeletionQV(Data::Orientation orientation, bool aligned,
                                          bool exciseSoftClips) const
{
    return FetchQualities(BamRecordTag::DELETION_QV, orientation, aligned, exciseSoftClips);
}

BamRecord& BamRecord::DeletionQV(const Data::QualityValues& deletionQVs)
{
    CreateOrEdit(BamRecordTag::DELETION_QV, deletionQVs.Fastq(), &impl_);
    return *this;
}

std::string BamRecord::DeletionTag(Data::Orientation orientation, bool aligned,
                                   bool exciseSoftClips) const
{
    return FetchBases(BamRecordTag::DELETION_TAG, orientation, aligned, exciseSoftClips);
}

BamRecord& BamRecord::DeletionTag(const std::string& tags)
{
    CreateOrEdit(BamRecordTag::DELETION_TAG, tags, &impl_);
    return *this;
}

std::vector<uint16_t> BamRecord::EncodePhotons(const std::vector<float>& data)
{
    std::vector<uint16_t> encoded;
    encoded.reserve(data.size());
    for (const auto& d : data) {
        encoded.emplace_back(d * photonFactor);
    }
    return encoded;
}

int BamRecord::EstimatedBytesUsed() const noexcept
{
    int result = impl_.EstimatedBytesUsed();
    result += sizeof(std::shared_ptr<BamHeader>);
    result += (2 * sizeof(Data::Position));
    result += sizeof(std::unique_ptr<Pulse2BaseCache>);
    if (p2bCache_) {
        result += p2bCache_->EstimatedBytesUsed();
    }
    return result;
}

std::string BamRecord::FetchBasesRaw(const BamRecordTag tag) const
{
    const Tag seqTag = impl_.TagValue(tag);
    return seqTag.ToString();
}

std::string BamRecord::FetchBases(const BamRecordTag tag, const Data::Orientation orientation,
                                  const bool aligned, const bool exciseSoftClips,
                                  const PulseBehavior pulseBehavior) const
{
    const bool isBamSeq = (tag == BamRecordTag::SEQ);
    const bool isPulse = BamRecordTags::IsPulse(tag);

    // fetch raw
    std::string bases;
    Data::Orientation current;
    if (isBamSeq) {  // SEQ stored in genomic orientation
        bases = impl_.Sequence();
        current = Data::Orientation::GENOMIC;
    } else {  // all tags stored in native orientation
        bases = FetchBasesRaw(tag);
        current = Data::Orientation::NATIVE;
    }

    // maybe strip 'squashed' pulse loci
    if (isPulse && pulseBehavior == PulseBehavior::BASECALLS_ONLY) {
        CalculatePulse2BaseCache();
        bases = p2bCache_->RemoveSquashedPulses(bases);
    }

    // if we need to touch CIGAR
    if (aligned || exciseSoftClips) {

        if (isPulse && pulseBehavior != PulseBehavior::BASECALLS_ONLY) {
            throw std::runtime_error{
                "[pbbam] BAM record ERROR: cannot return data at all pulses when gapping and/or "
                "soft-clipping are requested. Use PulseBehavior::BASECALLS_ONLY instead."};
        }

        // force into genomic orientation
        OrientBasesAsRequested(&bases, current, Data::Orientation::GENOMIC, impl_.IsReverseStrand(),
                               isPulse);
        current = Data::Orientation::GENOMIC;

        // clip & gapify as requested
        ClipAndGapifyBases(impl_, aligned, exciseSoftClips, &bases);
    }

    // return in the orientation requested
    OrientBasesAsRequested(&bases, current, orientation, impl_.IsReverseStrand(), isPulse);
    return bases;
}

Data::Frames BamRecord::FetchFramesRaw(const BamRecordTag tag) const
{
    const auto frameTag = impl_.TagValue(tag);
    if (frameTag.IsNull()) {
        return {};
    }  // throw ?

    // lossy frame codes
    if (frameTag.IsUInt8Array()) {
        const auto& decoder = [&]() -> Data::FrameEncoder {
            if (BamRecordTags::IsIPD(tag)) {
                return IpdEncoder(*this);
            } else {
                assert(BamRecordTags::IsPW(tag));
                return PwEncoder(*this);
            }
        }();
        return decoder.Decode(frameTag.ToUInt8Array());
    }

    // lossless frame data
    else {
        assert(frameTag.IsUInt16Array());
        return Data::Frames{frameTag.ToUInt16Array()};
    }
}

Data::Frames BamRecord::FetchFrames(const BamRecordTag tag, const Data::Orientation orientation,
                                    const bool aligned, const bool exciseSoftClips,
                                    const PulseBehavior pulseBehavior) const
{
    const bool isPulse = BamRecordTags::IsPulse(tag);

    // fetch raw
    Data::Frames frames = FetchFramesRaw(tag);
    Data::Orientation current = Data::Orientation::NATIVE;

    // maybe strip 'squashed' pulse loci
    if (isPulse && pulseBehavior == PulseBehavior::BASECALLS_ONLY) {
        CalculatePulse2BaseCache();
        frames.DataRaw() = p2bCache_->RemoveSquashedPulses(frames.Data());
    }

    // if we need to touch the CIGAR
    if (aligned || exciseSoftClips) {

        if (isPulse && pulseBehavior != PulseBehavior::BASECALLS_ONLY) {
            throw std::runtime_error{
                "[pbbam] BAM record ERROR: cannot return data at all pulses when gapping and/or "
                "soft-clipping are requested. Use PulseBehavior::BASECALLS_ONLY instead."};
        }

        // force into genomic orientation
        OrientTagDataAsRequested(&frames, current, Data::Orientation::GENOMIC,
                                 impl_.IsReverseStrand());
        current = Data::Orientation::GENOMIC;

        // clip & gapify as requested
        ClipAndGapifyFrames(impl_, aligned, exciseSoftClips, &frames);
    }

    // return in the orientation requested
    OrientTagDataAsRequested(&frames, current, orientation, impl_.IsReverseStrand());
    return frames;
}

std::vector<float> BamRecord::FetchPhotonsRaw(const BamRecordTag tag) const
{
    const auto frameTag = impl_.TagValue(tag);
    if (frameTag.IsNull()) {
        return {};
    }
    if (!frameTag.IsUInt16Array()) {
        throw std::runtime_error{
            "[pbbam] BAM record ERROR: photons are not a uint16_t array, tag " +
            BamRecordTags::LabelFor(tag)};
    }

    const auto data = frameTag.ToUInt16Array();
    std::vector<float> photons;
    photons.reserve(data.size());
    for (const auto& d : data) {
        photons.emplace_back(d / photonFactor);
    }
    return photons;
}

std::vector<float> BamRecord::FetchPhotons(const BamRecordTag tag,
                                           const Data::Orientation orientation, const bool aligned,
                                           const bool exciseSoftClips,
                                           const PulseBehavior pulseBehavior) const
{
    const bool isPulse = BamRecordTags::IsPulse(tag);

    // fetch raw
    auto data = FetchPhotonsRaw(tag);
    Data::Orientation current = Data::Orientation::NATIVE;

    if (isPulse && pulseBehavior == PulseBehavior::BASECALLS_ONLY) {
        // strip 'squashed' pulse loci
        CalculatePulse2BaseCache();
        data = p2bCache_->RemoveSquashedPulses(data);
    }

    if (aligned || exciseSoftClips) {

        if (isPulse && pulseBehavior != PulseBehavior::BASECALLS_ONLY) {
            throw std::runtime_error{
                "[pbbam] BAM record ERROR: cannot return data at all pulses when gapping and/or "
                "soft-clipping are requested. Use PulseBehavior::BASECALLS_ONLY instead."};
        }

        // force into genomic orientation
        OrientTagDataAsRequested(&data, current, Data::Orientation::GENOMIC,
                                 impl_.IsReverseStrand());
        current = Data::Orientation::GENOMIC;

        // clip & gapify as requested
        ClipAndGapifyPhotons(impl_, aligned, exciseSoftClips, &data);
    }

    // return in the orientation requested
    OrientTagDataAsRequested(&data, current, orientation, impl_.IsReverseStrand());
    return data;
}

Data::QualityValues BamRecord::FetchQualitiesRaw(const BamRecordTag tag) const
{
    const auto qvsTag = impl_.TagValue(tag);
    return Data::QualityValues::FromFastq(qvsTag.ToString());
}

Data::QualityValues BamRecord::FetchQualities(const BamRecordTag tag,
                                              const Data::Orientation orientation,
                                              const bool aligned, const bool exciseSoftClips,
                                              const PulseBehavior pulseBehavior) const
{
    // requested data info
    const bool isBamQual = (tag == BamRecordTag::QUAL);
    const bool isPulse = BamRecordTags::IsPulse(tag);

    // fetch raw
    Data::QualityValues quals;
    Data::Orientation current;
    if (isBamQual) {  // QUAL stored in genomic orientation
        quals = impl_.Qualities();
        current = Data::Orientation::GENOMIC;
    } else {  // all tags stored in native orientation
        quals = FetchQualitiesRaw(tag);
        current = Data::Orientation::NATIVE;
    }

    if (isPulse && pulseBehavior == PulseBehavior::BASECALLS_ONLY) {
        // strip 'squashed' pulse loci
        CalculatePulse2BaseCache();
        quals = p2bCache_->RemoveSquashedPulses(quals);
    }

    // if we need to touch CIGAR
    if (aligned || exciseSoftClips) {

        if (isPulse && pulseBehavior != PulseBehavior::BASECALLS_ONLY) {
            throw std::runtime_error{
                "[pbbam] BAM record ERROR: cannot return data at all pulses when gapping and/or "
                "soft-clipping are requested. Use PulseBehavior::BASECALLS_ONLY instead."};
        }

        // force into genomic orientation
        OrientTagDataAsRequested(&quals, current, Data::Orientation::GENOMIC,
                                 impl_.IsReverseStrand());
        current = Data::Orientation::GENOMIC;

        // clip & gapify as requested
        ClipAndGapifyQualities(impl_, aligned, exciseSoftClips, &quals);
    }

    // return in the orientation requested
    OrientTagDataAsRequested(&quals, current, orientation, impl_.IsReverseStrand());
    return quals;
}

std::vector<uint32_t> BamRecord::FetchUInt32sRaw(const BamRecordTag tag) const
{
    // fetch tag data
    const auto frameTag = impl_.TagValue(tag);
    if (frameTag.IsNull()) {
        return {};
    }
    if (!frameTag.IsUInt32Array()) {
        throw std::runtime_error{
            "[pbbam] BAM record ERROR: tag data are not a uint32_t array, tag " +
            BamRecordTags::LabelFor(tag)};
    }
    return frameTag.ToUInt32Array();
}

std::vector<uint32_t> BamRecord::FetchUInt32s(const BamRecordTag tag,
                                              const Data::Orientation orientation,
                                              const bool aligned, const bool exciseSoftClips,
                                              const PulseBehavior pulseBehavior) const
{
    const bool isPulse = BamRecordTags::IsPulse(tag);

    // fetch raw
    auto arr = FetchUInt32sRaw(tag);
    Data::Orientation current = Data::Orientation::NATIVE;

    if (isPulse && pulseBehavior == PulseBehavior::BASECALLS_ONLY) {
        // strip 'squashed' pulse loci
        CalculatePulse2BaseCache();
        arr = p2bCache_->RemoveSquashedPulses(arr);
    }

    if (aligned || exciseSoftClips) {

        if (isPulse && pulseBehavior != PulseBehavior::BASECALLS_ONLY) {
            throw std::runtime_error{
                "[pbbam] BAM record ERROR: cannot return data at all pulses when gapping and/or "
                "soft-clipping are requested. Use PulseBehavior::BASECALLS_ONLY instead."};
        }

        // force into genomic orientation
        OrientTagDataAsRequested(&arr, current, Data::Orientation::GENOMIC,
                                 impl_.IsReverseStrand());
        current = Data::Orientation::GENOMIC;

        // clip & gapify as requested
        ClipAndGapifyUInts(impl_, aligned, exciseSoftClips, &arr);
    }

    // return in the orientation requested
    OrientTagDataAsRequested(&arr, current, orientation, impl_.IsReverseStrand());
    return arr;
}

std::vector<uint8_t> BamRecord::FetchUInt8sRaw(const BamRecordTag tag) const
{
    // fetch tag data
    const auto frameTag = impl_.TagValue(tag);
    if (frameTag.IsNull()) {
        return {};
    }
    if (!frameTag.IsUInt8Array()) {
        throw std::runtime_error{
            "[pbbam] BAM record ERROR: tag data are not a uint8_t array, tag " +
            BamRecordTags::LabelFor(tag)};
    }
    return frameTag.ToUInt8Array();
}

std::vector<uint8_t> BamRecord::FetchUInt8s(const BamRecordTag tag,
                                            const Data::Orientation orientation, const bool aligned,
                                            const bool exciseSoftClips,
                                            const PulseBehavior pulseBehavior) const
{
    const bool isPulse = BamRecordTags::IsPulse(tag);

    // fetch raw
    auto arr = FetchUInt8sRaw(tag);
    Data::Orientation current = Data::Orientation::NATIVE;

    if (isPulse && pulseBehavior == PulseBehavior::BASECALLS_ONLY) {
        // strip 'squashed' pulse loci
        CalculatePulse2BaseCache();
        arr = p2bCache_->RemoveSquashedPulses(arr);
    }

    if (aligned || exciseSoftClips) {

        if (isPulse && pulseBehavior != PulseBehavior::BASECALLS_ONLY) {
            throw std::runtime_error{
                "[pbbam] BAM record ERROR: cannot return data at all pulses when gapping and/or "
                "soft-clipping are requested. Use PulseBehavior::BASECALLS_ONLY instead."};
        }

        // force into genomic orientation
        OrientTagDataAsRequested(&arr, current, Data::Orientation::GENOMIC,
                                 impl_.IsReverseStrand());
        current = Data::Orientation::GENOMIC;

        // clip & gapify as requested
        ClipAndGapifyUInt8s(impl_, aligned, exciseSoftClips, &arr);
    }

    // return in the orientation requested
    OrientTagDataAsRequested(&arr, current, orientation, impl_.IsReverseStrand());
    return arr;
}

Data::Frames BamRecord::ForwardIPD(Data::Orientation orientation, bool aligned,
                                   bool exciseSoftClips) const
{
    return FetchFrames(BamRecordTag::FORWARD_IPD, orientation, aligned, exciseSoftClips);
}

BamRecord& BamRecord::ForwardIPD(const Data::Frames& frames, const Data::FrameCodec encoding)
{
    const auto& frameData = frames.Data();
    if (encoding == Data::FrameCodec::RAW) {
        CreateOrEdit(BamRecordTag::FORWARD_IPD, frameData, &impl_);
    } else {
        const auto encoder = IpdEncoder(*this);
        CreateOrEdit(BamRecordTag::FORWARD_IPD, encoder.Encode(frameData), &impl_);
    }
    return *this;
}

Data::Frames BamRecord::ForwardPulseWidth(Data::Orientation orientation, bool aligned,
                                          bool exciseSoftClips) const
{
    return FetchFrames(BamRecordTag::FORWARD_PW, orientation, aligned, exciseSoftClips);
}

BamRecord& BamRecord::ForwardPulseWidth(const Data::Frames& frames, const Data::FrameCodec encoding)
{
    const auto& frameData = frames.Data();
    if (encoding == Data::FrameCodec::RAW) {
        CreateOrEdit(BamRecordTag::FORWARD_PW, frameData, &impl_);
    } else {
        const auto encoder = PwEncoder(*this);
        CreateOrEdit(BamRecordTag::FORWARD_PW, encoder.Encode(frameData), &impl_);
    }
    return *this;
}

std::string BamRecord::FullName() const { return impl_.Name(); }

bool BamRecord::HasAltLabelQV() const { return impl_.HasTag(BamRecordTag::ALT_LABEL_QV); }

bool BamRecord::HasAltLabelTag() const { return impl_.HasTag(BamRecordTag::ALT_LABEL_TAG); }

bool BamRecord::HasBarcodes() const { return impl_.HasTag(BamRecordTag::BARCODES); }

bool BamRecord::HasBarcodeQuality() const { return impl_.HasTag(BamRecordTag::BARCODE_QUALITY); }

bool BamRecord::HasLabelQV() const { return impl_.HasTag(BamRecordTag::LABEL_QV); }

bool BamRecord::HasDeletionQV() const { return impl_.HasTag(BamRecordTag::DELETION_QV); }

bool BamRecord::HasDeletionTag() const { return impl_.HasTag(BamRecordTag::DELETION_TAG); }

bool BamRecord::HasForwardIPD() const
{
    return impl_.HasTag(BamRecordTag::FORWARD_IPD) &&
           !impl_.TagValue(BamRecordTag::FORWARD_IPD).IsNull();
}

bool BamRecord::HasForwardPulseWidth() const
{
    return impl_.HasTag(BamRecordTag::FORWARD_PW) &&
           !impl_.TagValue(BamRecordTag::FORWARD_PW).IsNull();
}

bool BamRecord::HasHoleNumber() const
{
    return impl_.HasTag(BamRecordTag::HOLE_NUMBER) &&
           !impl_.TagValue(BamRecordTag::HOLE_NUMBER).IsNull();
}

bool BamRecord::HasInsertionQV() const { return impl_.HasTag(BamRecordTag::INSERTION_QV); }

bool BamRecord::HasNumPasses() const { return impl_.HasTag(BamRecordTag::NUM_PASSES); }

bool BamRecord::HasPreBaseFrames() const { return HasIPD(); }

bool BamRecord::HasIPD() const { return impl_.HasTag(BamRecordTag::IPD); }

bool BamRecord::HasLocalContextFlags() const { return impl_.HasTag(BamRecordTag::CONTEXT_FLAGS); }

bool BamRecord::HasMergeQV() const { return impl_.HasTag(BamRecordTag::MERGE_QV); }

bool BamRecord::HasPulseMergeQV() const { return impl_.HasTag(BamRecordTag::PULSE_MERGE_QV); }

bool BamRecord::HasPkmean() const { return impl_.HasTag(BamRecordTag::PKMEAN); }

bool BamRecord::HasPkmean2() const { return impl_.HasTag(BamRecordTag::PKMEAN_2); }

bool BamRecord::HasPkmid() const { return impl_.HasTag(BamRecordTag::PKMID); }

bool BamRecord::HasPkmid2() const { return impl_.HasTag(BamRecordTag::PKMID_2); }

bool BamRecord::HasPrePulseFrames() const { return impl_.HasTag(BamRecordTag::PRE_PULSE_FRAMES); }

bool BamRecord::HasPulseCall() const
{
    return impl_.HasTag(BamRecordTag::PULSE_CALL) &&
           !impl_.TagValue(BamRecordTag::PULSE_CALL).IsNull();
}

bool BamRecord::HasPulseExclusion() const { return impl_.HasTag(BamRecordTag::PULSE_EXCLUSION); }

bool BamRecord::HasPulseCallWidth() const { return impl_.HasTag(BamRecordTag::PULSE_CALL_WIDTH); }

bool BamRecord::HasPulseWidth() const { return impl_.HasTag(BamRecordTag::PULSE_WIDTH); }

bool BamRecord::HasQueryEnd() const { return impl_.HasTag(BamRecordTag::QUERY_END); }

bool BamRecord::HasQueryEndFrameNumber() const
{
    return impl_.HasTag(BamRecordTag::QUERY_END_FRAME_NUMBER);
}

bool BamRecord::HasQueryStart() const { return impl_.HasTag(BamRecordTag::QUERY_START); }

bool BamRecord::HasQueryStartFrameNumber() const
{
    return impl_.HasTag(BamRecordTag::QUERY_START_FRAME_NUMBER);
}

bool BamRecord::HasReadAccuracy() const
{
    return impl_.HasTag(BamRecordTag::READ_ACCURACY) &&
           !impl_.TagValue(BamRecordTag::READ_ACCURACY).IsNull();
}

bool BamRecord::HasReverseIPD() const
{
    return impl_.HasTag(BamRecordTag::REVERSE_IPD) &&
           !impl_.TagValue(BamRecordTag::REVERSE_IPD).IsNull();
}

bool BamRecord::HasReversePulseWidth() const
{
    return impl_.HasTag(BamRecordTag::REVERSE_PW) &&
           !impl_.TagValue(BamRecordTag::REVERSE_PW).IsNull();
}

bool BamRecord::HasScrapRegionType() const
{
    return impl_.HasTag(BamRecordTag::SCRAP_REGION_TYPE) &&
           !impl_.TagValue(BamRecordTag::SCRAP_REGION_TYPE).IsNull();
}

bool BamRecord::HasScrapZmwType() const
{
    return impl_.HasTag(BamRecordTag::SCRAP_ZMW_TYPE) &&
           !impl_.TagValue(BamRecordTag::SCRAP_ZMW_TYPE).IsNull();
}

bool BamRecord::HasStartFrame() const { return impl_.HasTag(BamRecordTag::START_FRAME); }

bool BamRecord::HasSignalToNoise() const { return impl_.HasTag(BamRecordTag::SIGNAL_TO_NOISE); }

bool BamRecord::HasSubstitutionQV() const { return impl_.HasTag(BamRecordTag::SUBSTITUTION_QV); }

bool BamRecord::HasSubstitutionTag() const { return impl_.HasTag(BamRecordTag::SUBSTITUTION_TAG); }

BamHeader BamRecord::Header() const { return header_; }

int32_t BamRecord::HoleNumber() const
{
    const Tag holeNumber = impl_.TagValue(BamRecordTag::HOLE_NUMBER);
    if (!holeNumber.IsNull()) {
        return holeNumber.ToInt32();
    }

    // missing zm tag - try to pull from name
    return HoleNumberFromName(FullName());
}

BamRecord& BamRecord::HoleNumber(const int32_t holeNumber)
{
    CreateOrEdit(BamRecordTag::HOLE_NUMBER, holeNumber, &impl_);
    return *this;
}

BamRecordImpl& BamRecord::Impl() { return impl_; }

const BamRecordImpl& BamRecord::Impl() const { return impl_; }

Data::QualityValues BamRecord::InsertionQV(Data::Orientation orientation, bool aligned,
                                           bool exciseSoftClips) const
{
    return FetchQualities(BamRecordTag::INSERTION_QV, orientation, aligned, exciseSoftClips);
}

BamRecord& BamRecord::InsertionQV(const Data::QualityValues& insertionQVs)
{
    CreateOrEdit(BamRecordTag::INSERTION_QV, insertionQVs.Fastq(), &impl_);
    return *this;
}

Data::Frames BamRecord::IPD(Data::Orientation orientation, bool aligned, bool exciseSoftClips) const
{
    return FetchFrames(BamRecordTag::IPD, orientation, aligned, exciseSoftClips);
}

BamRecord& BamRecord::IPD(const Data::Frames& frames, const FrameEncodingType encoding)
{
    if (encoding == FrameEncodingType::LOSSY) {
        return IPD(frames, ReadGroup().IpdCodec());
    } else {
        return IPD(frames, Data::FrameCodec::RAW);
    }
}

BamRecord& BamRecord::IPD(const Data::Frames& frames, const Data::FrameCodec encoding)
{
    const auto& frameData = frames.Data();
    if (encoding == Data::FrameCodec::RAW) {
        CreateOrEdit(BamRecordTag::IPD, frameData, &impl_);
    } else {
        const auto encoder = IpdEncoder(*this);
        CreateOrEdit(BamRecordTag::IPD, encoder.Encode(frameData), &impl_);
    }
    return *this;
}

Data::Frames BamRecord::IPDRaw(Data::Orientation orientation) const
{
    const auto tagName = BamRecordTags::LabelFor(BamRecordTag::IPD);
    const auto frameTag = impl_.TagValue(tagName);
    if (frameTag.IsNull()) {
        return {};
    }

    Data::Frames frames;

    // lossy frame codes
    if (frameTag.IsUInt8Array()) {
        const auto codes = frameTag.ToUInt8Array();
        const std::vector<uint16_t> codes16(codes.begin(), codes.end());
        frames.Data(std::move(codes16));
    }

    // lossless frame data
    else {
        assert(frameTag.IsUInt16Array());
        frames.Data(frameTag.ToUInt16Array());
    }

    // return in requested orientation
    OrientTagDataAsRequested(&frames,
                             Data::Orientation::NATIVE,  // current
                             orientation,                // requested
                             impl_.IsReverseStrand());
    return frames;
}

bool BamRecord::IsMapped() const { return impl_.IsMapped(); }

Data::QualityValues BamRecord::LabelQV(Data::Orientation orientation, bool aligned,
                                       bool exciseSoftClips, PulseBehavior pulseBehavior) const
{
    return FetchQualities(BamRecordTag::LABEL_QV, orientation, aligned, exciseSoftClips,
                          pulseBehavior);
}

BamRecord& BamRecord::LabelQV(const Data::QualityValues& labelQVs)
{
    CreateOrEdit(BamRecordTag::LABEL_QV, labelQVs.Fastq(), &impl_);
    return *this;
}

Data::LocalContextFlags BamRecord::LocalContextFlags() const
{
    assert(HasLocalContextFlags());
    const auto tagName = BamRecordTags::LabelFor(BamRecordTag::CONTEXT_FLAGS);
    const auto cxTag = impl_.TagValue(tagName);
    return static_cast<Data::LocalContextFlags>(cxTag.ToUInt8());
}

BamRecord& BamRecord::LocalContextFlags(const Data::LocalContextFlags flags)
{
    CreateOrEdit(BamRecordTag::CONTEXT_FLAGS, static_cast<uint8_t>(flags), &impl_);
    return *this;
}

BamRecord& BamRecord::Map(const int32_t referenceId, const Data::Position refStart,
                          const Data::Strand strand, const Data::Cigar& cigar,
                          const uint8_t mappingQuality)
{
    impl_.Position(refStart);
    impl_.ReferenceId(referenceId);
    impl_.CigarData(cigar);
    impl_.MapQuality(mappingQuality);
    impl_.SetMapped(true);

    if (strand == Data::Strand::FORWARD) {
        impl_.SetReverseStrand(false);

    } else {
        assert(strand == Data::Strand::REVERSE);
        impl_.SetReverseStrand(true);

        // switch seq & qual
        auto sequence = impl_.Sequence();
        auto qualities = impl_.Qualities();

        ReverseComplement(sequence);
        Reverse(qualities);

        impl_.SetSequenceAndQualities(sequence, qualities.Fastq());
    }

    // reset any cached aligned start/end
    alignedStart_ = Data::UnmappedPosition;
    alignedEnd_ = Data::UnmappedPosition;

    return *this;
}

BamRecord BamRecord::Mapped(const BamRecord& input, const int32_t referenceId,
                            const Data::Position refStart, const Data::Strand strand,
                            const Data::Cigar& cigar, const uint8_t mappingQuality)
{
    return input.Mapped(referenceId, refStart, strand, cigar, mappingQuality);
}

BamRecord BamRecord::Mapped(const int32_t referenceId, const Data::Position refStart,
                            const Data::Strand strand, const Data::Cigar& cigar,
                            const uint8_t mappingQuality) const
{
    BamRecord result(*this);
    result.Map(referenceId, refStart, strand, cigar, mappingQuality);
    return result;
}

uint8_t BamRecord::MapQuality() const { return impl_.MapQuality(); }

Data::QualityValues BamRecord::MergeQV(Data::Orientation orientation, bool aligned,
                                       bool exciseSoftClips) const
{
    return FetchQualities(BamRecordTag::MERGE_QV, orientation, aligned, exciseSoftClips);
}

BamRecord& BamRecord::MergeQV(const Data::QualityValues& mergeQVs)
{
    CreateOrEdit(BamRecordTag::MERGE_QV, mergeQVs.Fastq(), &impl_);
    return *this;
}

std::string BamRecord::MovieName() const
{
    const auto& rgId = ReadGroupId();
    if (!rgId.empty()) {
        return header_.ReadGroup(rgId).MovieName();
    } else {
        const auto nameParts = Split(FullName(), '/');
        if (nameParts.empty()) {
            throw std::runtime_error{"[pbbam] BAM record ERROR: has malformed name: '" +
                                     FullName() + "'"};
        }
        return nameParts[0];
    }
}

size_t BamRecord::NumDeletedBases() const { return NumInsertedAndDeletedBases().second; }

size_t BamRecord::NumDeletionOperations() const
{
    return NumInsertionAndDeletionOperations().second;
}

std::pair<size_t, size_t> BamRecord::NumInsertedAndDeletedBases() const
{
    size_t nInsBases = 0;
    size_t nDelBases = 0;

    auto& b = BamRecordMemory::GetRawData(this);
    uint32_t* cigarData = bam_get_cigar(b.get());
    for (uint32_t i = 0; i < b->core.n_cigar; ++i) {
        const auto type = static_cast<Data::CigarOperationType>(bam_cigar_op(cigarData[i]));
        if (type == Data::CigarOperationType::INSERTION) {
            nInsBases += bam_cigar_oplen(cigarData[i]);
        } else if (type == Data::CigarOperationType::DELETION) {
            nDelBases += bam_cigar_oplen(cigarData[i]);
        }
    }
    return {nInsBases, nDelBases};
}

size_t BamRecord::NumInsertedBases() const { return NumInsertedAndDeletedBases().first; }

std::pair<size_t, size_t> BamRecord::NumInsertionAndDeletionOperations() const
{
    size_t nInsOps = 0;
    size_t nDelOps = 0;

    auto& b = BamRecordMemory::GetRawData(this);
    uint32_t* cigarData = bam_get_cigar(b.get());
    for (uint32_t i = 0; i < b->core.n_cigar; ++i) {
        const auto type = static_cast<Data::CigarOperationType>(bam_cigar_op(cigarData[i]));
        if (type == Data::CigarOperationType::INSERTION) {
            ++nInsOps;
        } else if (type == Data::CigarOperationType::DELETION) {
            ++nDelOps;
        }
    }
    return {nInsOps, nDelOps};
}

size_t BamRecord::NumInsertionOperations() const
{
    return NumInsertionAndDeletionOperations().first;
}

size_t BamRecord::NumMatches() const { return NumMatchesAndMismatches().first; }

std::pair<size_t, size_t> BamRecord::NumMatchesAndMismatches() const
{
    std::pair<size_t, size_t> result = std::make_pair(0, 0);

    auto& b = BamRecordMemory::GetRawData(this);
    uint32_t* cigarData = bam_get_cigar(b.get());
    for (uint32_t i = 0; i < b->core.n_cigar; ++i) {
        const auto type = static_cast<Data::CigarOperationType>(bam_cigar_op(cigarData[i]));
        if (type == Data::CigarOperationType::SEQUENCE_MATCH) {
            result.first += bam_cigar_oplen(cigarData[i]);
        } else if (type == Data::CigarOperationType::SEQUENCE_MISMATCH) {
            result.second += bam_cigar_oplen(cigarData[i]);
        }
    }
    return result;
}

size_t BamRecord::NumMismatches() const { return NumMatchesAndMismatches().second; }

int32_t BamRecord::NumPasses() const
{
    const auto tagName = BamRecordTags::LabelFor(BamRecordTag::NUM_PASSES);
    const auto numPasses = impl_.TagValue(tagName);
    return numPasses.ToInt32();
}

BamRecord& BamRecord::NumPasses(const int32_t numPasses)
{
    CreateOrEdit(BamRecordTag::NUM_PASSES, numPasses, &impl_);
    return *this;
}

std::vector<float> BamRecord::Pkmean(Data::Orientation orientation, bool aligned,
                                     bool exciseSoftClips, PulseBehavior pulseBehavior) const
{
    return FetchPhotons(BamRecordTag::PKMEAN, orientation, aligned, exciseSoftClips, pulseBehavior);
}

BamRecord& BamRecord::Pkmean(const std::vector<float>& photons)
{
    Pkmean(EncodePhotons(photons));
    return *this;
}

BamRecord& BamRecord::Pkmean(const std::vector<uint16_t>& encodedPhotons)
{
    CreateOrEdit(BamRecordTag::PKMEAN, encodedPhotons, &impl_);
    return *this;
}

std::vector<float> BamRecord::Pkmid(Data::Orientation orientation, bool aligned,
                                    bool exciseSoftClips, PulseBehavior pulseBehavior) const
{
    return FetchPhotons(BamRecordTag::PKMID, orientation, aligned, exciseSoftClips, pulseBehavior);
}

BamRecord& BamRecord::Pkmid(const std::vector<float>& photons)
{
    Pkmid(EncodePhotons(photons));
    return *this;
}

BamRecord& BamRecord::Pkmid(const std::vector<uint16_t>& encodedPhotons)
{
    CreateOrEdit(BamRecordTag::PKMID, encodedPhotons, &impl_);
    return *this;
}

std::vector<float> BamRecord::Pkmean2(Data::Orientation orientation, bool aligned,
                                      bool exciseSoftClips, PulseBehavior pulseBehavior) const
{
    return FetchPhotons(BamRecordTag::PKMEAN_2, orientation, aligned, exciseSoftClips,
                        pulseBehavior);
}

BamRecord& BamRecord::Pkmean2(const std::vector<float>& photons)
{
    Pkmean2(EncodePhotons(photons));
    return *this;
}

BamRecord& BamRecord::Pkmean2(const std::vector<uint16_t>& encodedPhotons)
{
    CreateOrEdit(BamRecordTag::PKMEAN_2, encodedPhotons, &impl_);
    return *this;
}

std::vector<float> BamRecord::Pkmid2(Data::Orientation orientation, bool aligned,
                                     bool exciseSoftClips, PulseBehavior pulseBehavior) const
{
    return FetchPhotons(BamRecordTag::PKMID_2, orientation, aligned, exciseSoftClips,
                        pulseBehavior);
}

BamRecord& BamRecord::Pkmid2(const std::vector<float>& photons)
{
    Pkmid2(EncodePhotons(photons));
    return *this;
}

BamRecord& BamRecord::Pkmid2(const std::vector<uint16_t>& encodedPhotons)
{
    CreateOrEdit(BamRecordTag::PKMID_2, encodedPhotons, &impl_);
    return *this;
}

Data::Frames BamRecord::PreBaseFrames(Data::Orientation orientation, bool aligned,
                                      bool exciseSoftClips) const
{
    return IPD(orientation, aligned, exciseSoftClips);
}

BamRecord& BamRecord::PreBaseFrames(const Data::Frames& frames, const FrameEncodingType encoding)
{
    return IPD(frames, encoding);
}

BamRecord& BamRecord::PreBaseFrames(const Data::Frames& frames, const Data::FrameCodec encoding)
{
    return IPD(frames, encoding);
}

Data::Frames BamRecord::PrePulseFrames(Data::Orientation orientation, bool aligned,
                                       bool exciseSoftClips, PulseBehavior pulseBehavior) const
{
    return FetchFrames(BamRecordTag::PRE_PULSE_FRAMES, orientation, aligned, exciseSoftClips,
                       pulseBehavior);
}

BamRecord& BamRecord::PrePulseFrames(const Data::Frames& frames, const FrameEncodingType encoding)
{
    if (encoding == FrameEncodingType::LOSSY) {
        return PrePulseFrames(frames, Data::FrameCodec::V1);
    } else {
        return PrePulseFrames(frames, Data::FrameCodec::RAW);
    }
}

BamRecord& BamRecord::PrePulseFrames(const Data::Frames& frames, const Data::FrameCodec encoding)
{
    if (encoding == Data::FrameCodec::V1) {
        CreateOrEdit(BamRecordTag::PRE_PULSE_FRAMES, frames.Encode(), &impl_);
    } else {
        CreateOrEdit(BamRecordTag::PRE_PULSE_FRAMES, frames.Data(), &impl_);
    }
    return *this;
}

Data::Frames BamRecord::PulseWidthRaw(Data::Orientation orientation, bool /* aligned */,
                                      bool /* exciseSoftClips */) const
{
    const auto tagName = BamRecordTags::LabelFor(BamRecordTag::PULSE_WIDTH);
    const auto frameTag = impl_.TagValue(tagName);
    if (frameTag.IsNull()) {
        return {};
    }

    Data::Frames frames;

    // lossy frame codes
    if (frameTag.IsUInt8Array()) {
        const auto codes = frameTag.ToUInt8Array();
        const std::vector<uint16_t> codes16(codes.begin(), codes.end());
        frames.Data(std::move(codes16));
    }

    // lossless frame data
    else {
        assert(frameTag.IsUInt16Array());
        frames.Data(frameTag.ToUInt16Array());
    }

    // return in requested orientation
    OrientTagDataAsRequested(&frames,
                             Data::Orientation::NATIVE,  // current
                             orientation,                // requested
                             impl_.IsReverseStrand());
    return frames;
}

Data::QualityValues BamRecord::PulseMergeQV(Data::Orientation orientation, bool aligned,
                                            bool exciseSoftClips, PulseBehavior pulseBehavior) const
{
    return FetchQualities(BamRecordTag::PULSE_MERGE_QV, orientation, aligned, exciseSoftClips,
                          pulseBehavior);
}

BamRecord& BamRecord::PulseMergeQV(const Data::QualityValues& mergeQVs)
{
    CreateOrEdit(BamRecordTag::PULSE_MERGE_QV, mergeQVs.Fastq(), &impl_);
    return *this;
}

std::string BamRecord::PulseCall(Data::Orientation orientation, bool aligned, bool exciseSoftClips,
                                 PulseBehavior pulseBehavior) const
{
    return FetchBases(BamRecordTag::PULSE_CALL, orientation, aligned, exciseSoftClips,
                      pulseBehavior);
}

BamRecord& BamRecord::PulseCall(const std::string& tags)
{
    CreateOrEdit(BamRecordTag::PULSE_CALL, tags, &impl_);
    return *this;
}

Data::Frames BamRecord::PulseCallWidth(Data::Orientation orientation, bool aligned,
                                       bool exciseSoftClips, PulseBehavior pulseBehavior) const
{
    return FetchFrames(BamRecordTag::PULSE_CALL_WIDTH, orientation, aligned, exciseSoftClips,
                       pulseBehavior);
}

BamRecord& BamRecord::PulseCallWidth(const Data::Frames& frames, const FrameEncodingType encoding)
{
    if (encoding == FrameEncodingType::LOSSY) {
        return PulseCallWidth(frames, Data::FrameCodec::V1);
    } else {
        return PulseCallWidth(frames, Data::FrameCodec::RAW);
    }
}

BamRecord& BamRecord::PulseCallWidth(const Data::Frames& frames, const Data::FrameCodec encoding)
{
    if (encoding == Data::FrameCodec::V1) {
        CreateOrEdit(BamRecordTag::PULSE_CALL_WIDTH, frames.Encode(), &impl_);
    } else {
        CreateOrEdit(BamRecordTag::PULSE_CALL_WIDTH, frames.Data(), &impl_);
    }
    return *this;
}

std::vector<BAM::PulseExclusionReason> BamRecord::PulseExclusionReason(
    Data::Orientation orientation, bool aligned, bool exciseSoftClips,
    PulseBehavior pulseBehavior) const
{
    std::vector<BAM::PulseExclusionReason> reasons;

    const auto reasonNums = FetchUInt8s(BamRecordTag::PULSE_EXCLUSION, orientation, aligned,
                                        exciseSoftClips, pulseBehavior);

    std::transform(reasonNums.cbegin(), reasonNums.cend(), std::back_inserter(reasons),
                   [](const uint8_t num) { return static_cast<BAM::PulseExclusionReason>(num); });

    return reasons;
}

BamRecord& BamRecord::PulseExclusionReason(const std::vector<BAM::PulseExclusionReason>& reasons)
{
    std::vector<uint8_t> reasonNums;
    std::transform(
        reasons.cbegin(), reasons.cend(), std::back_inserter(reasonNums),
        [](const BAM::PulseExclusionReason& reason) { return static_cast<uint8_t>(reason); });

    CreateOrEdit(BamRecordTag::PULSE_EXCLUSION, reasonNums, &impl_);
    return *this;
}

Data::Frames BamRecord::PulseWidth(Data::Orientation orientation, bool aligned,
                                   bool exciseSoftClips) const
{
    return FetchFrames(BamRecordTag::PULSE_WIDTH, orientation, aligned, exciseSoftClips,
                       PulseBehavior::ALL);
}

BamRecord& BamRecord::PulseWidth(const Data::Frames& frames, const FrameEncodingType encoding)
{
    if (encoding == FrameEncodingType::LOSSY) {
        return PulseWidth(frames, ReadGroup().PulseWidthCodec());
    } else {
        return PulseWidth(frames, Data::FrameCodec::RAW);
    }
}

BamRecord& BamRecord::PulseWidth(const Data::Frames& frames, const Data::FrameCodec encoding)
{
    const auto& frameData = frames.Data();
    if (encoding == Data::FrameCodec::RAW) {
        CreateOrEdit(BamRecordTag::PULSE_WIDTH, frameData, &impl_);
    } else {
        const auto encoder = PwEncoder(*this);
        CreateOrEdit(BamRecordTag::PULSE_WIDTH, encoder.Encode(frameData), &impl_);
    }
    return *this;
}

Data::QualityValues BamRecord::Qualities(Data::Orientation orientation, bool aligned,
                                         bool exciseSoftClips) const
{
    return FetchQualities(BamRecordTag::QUAL, orientation, aligned, exciseSoftClips);
}

Data::Position BamRecord::QueryEnd() const
{
    // try 'qe' tag
    const auto tagName = BamRecordTags::LabelFor(BamRecordTag::QUERY_END);
    const auto qe = impl_.TagValue(tagName);
    if (!qe.IsNull()) {
        return qe.ToInt32();
    }

    // tag missing, need to check movie name (fallback for non-PB BAMs, but ignore for CCS reads)
    RecordType type;
    try {
        type = Type();
    } catch (std::exception&) {
        return 0;
    }
    if (type == RecordType::CCS) {
        throw std::runtime_error{
            "[pbbam] BAM record ERROR: no query end is available for CCS read type"};
    }
    if (type == RecordType::TRANSCRIPT) {
        throw std::runtime_error{
            "[pbbam] BAM record ERROR: no query end is available for transcript read type"};
    }

    // PacBio BAM, non-CCS/transcript
    try {
        return QueryEndFromName(FullName());
    } catch (std::exception&) {
        // return fallback position
        return 0;
    }
}

BamRecord& BamRecord::QueryEnd(const Data::Position pos)
{
    CreateOrEdit(BamRecordTag::QUERY_END, static_cast<int32_t>(pos), &impl_);
    UpdateName();
    return *this;
}

int32_t BamRecord::QueryEndFrameNumber() const
{
    const auto tagName = BamRecordTags::LabelFor(BamRecordTag::QUERY_END_FRAME_NUMBER);
    const auto qs = impl_.TagValue(tagName);
    if (!qs.IsNull()) {
        return qs.ToInt32();
    }
    return 0;
}

BamRecord& BamRecord::QueryEndFrameNumber(const int32_t frameNumber)
{
    CreateOrEdit(BamRecordTag::QUERY_END_FRAME_NUMBER, frameNumber, &impl_);
    return *this;
}

Data::Position BamRecord::QueryStart() const
{
    // try 'qs' tag
    const auto tagName = BamRecordTags::LabelFor(BamRecordTag::QUERY_START);
    const auto qs = impl_.TagValue(tagName);
    if (!qs.IsNull()) {
        return qs.ToInt32();
    }

    // tag missing, need to check movie name (fallback for non-PB BAMs, but ignore for CCS reads)
    RecordType type;
    try {
        type = Type();
    } catch (std::exception&) {
        return 0;
    }
    if (type == RecordType::CCS) {
        throw std::runtime_error{
            "[pbbam] BAM record ERROR: no query start is available for CCS read type"};
    }
    if (type == RecordType::TRANSCRIPT) {
        throw std::runtime_error{
            "[pbbam] BAM record ERROR: no query start is available for transcript read type"};
    }

    // PacBio BAM, non-CCS/transcript
    try {
        return QueryStartFromName(FullName());
    } catch (std::exception&) {
        // return fallback position
        return 0;
    }
}

BamRecord& BamRecord::QueryStart(const Data::Position pos)
{
    CreateOrEdit(BamRecordTag::QUERY_START, static_cast<int32_t>(pos), &impl_);
    UpdateName();
    return *this;
}

int32_t BamRecord::QueryStartFrameNumber() const
{
    const auto tagName = BamRecordTags::LabelFor(BamRecordTag::QUERY_START_FRAME_NUMBER);
    const auto qs = impl_.TagValue(tagName);
    if (!qs.IsNull()) {
        return qs.ToInt32();
    }
    return 0;
}

BamRecord& BamRecord::QueryStartFrameNumber(const int32_t frameNumber)
{
    CreateOrEdit(BamRecordTag::QUERY_START_FRAME_NUMBER, frameNumber, &impl_);
    return *this;
}

Data::Accuracy BamRecord::ReadAccuracy() const
{
    const auto tagName = BamRecordTags::LabelFor(BamRecordTag::READ_ACCURACY);
    const auto readAccuracy = impl_.TagValue(tagName);
    return {readAccuracy.ToFloat()};
}

BamRecord& BamRecord::ReadAccuracy(const Data::Accuracy& accuracy)
{
    CreateOrEdit(BamRecordTag::READ_ACCURACY, static_cast<float>(accuracy), &impl_);
    return *this;
}

ReadGroupInfo BamRecord::ReadGroup() const { return header_.ReadGroup(ReadGroupId()); }

BamRecord& BamRecord::ReadGroup(const ReadGroupInfo& rg)
{
    CreateOrEdit(BamRecordTag::READ_GROUP, rg.Id(), &impl_);
    UpdateName();
    return *this;
}

std::string BamRecord::ReadGroupId() const
{
    const auto tagName = BamRecordTags::LabelFor(BamRecordTag::READ_GROUP);
    const auto rgTag = impl_.TagValue(tagName);
    if (rgTag.IsNull()) {
        return {};
    }
    return rgTag.ToString();
}

std::string BamRecord::ReadGroupBaseId() const { return ReadGroup().BaseId(); }

BamRecord& BamRecord::ReadGroupId(const std::string& id)
{
    CreateOrEdit(BamRecordTag::READ_GROUP, id, &impl_);
    UpdateName();
    return *this;
}

int32_t BamRecord::ReadGroupNumericId() const { return ReadGroupInfo::IdToInt(ReadGroupBaseId()); }

Data::Position BamRecord::ReferenceEnd() const
{
    if (!impl_.IsMapped()) {
        return Data::UnmappedPosition;
    }
    const auto& htsData = BamRecordMemory::GetRawData(impl_);
    if (!htsData) {
        return Data::UnmappedPosition;
    }
    return bam_endpos(htsData.get());
}

int32_t BamRecord::ReferenceId() const { return impl_.ReferenceId(); }

std::string BamRecord::ReferenceName() const
{
    if (IsMapped()) {
        return Header().SequenceName(ReferenceId());
    } else {
        throw std::runtime_error{
            "[pbbam] BAM record ERROR: unmapped record has no associated reference name"};
    }
}

Data::Position BamRecord::ReferenceStart() const { return impl_.Position(); }

Data::Frames BamRecord::ReverseIPD(Data::Orientation orientation, bool aligned,
                                   bool exciseSoftClips) const
{
    return FetchFrames(BamRecordTag::REVERSE_IPD, orientation, aligned, exciseSoftClips);
}

BamRecord& BamRecord::ReverseIPD(const Data::Frames& frames, const Data::FrameCodec encoding)
{
    const auto& frameData = frames.Data();
    if (encoding == Data::FrameCodec::RAW) {
        CreateOrEdit(BamRecordTag::REVERSE_IPD, frameData, &impl_);
    } else {
        const auto encoder = IpdEncoder(*this);
        CreateOrEdit(BamRecordTag::REVERSE_IPD, encoder.Encode(frameData), &impl_);
    }
    return *this;
}

Data::Frames BamRecord::ReversePulseWidth(Data::Orientation orientation, bool aligned,
                                          bool exciseSoftClips) const
{
    return FetchFrames(BamRecordTag::REVERSE_PW, orientation, aligned, exciseSoftClips);
}

BamRecord& BamRecord::ReversePulseWidth(const Data::Frames& frames, const Data::FrameCodec encoding)
{
    const auto& frameData = frames.Data();
    if (encoding == Data::FrameCodec::RAW) {
        CreateOrEdit(BamRecordTag::REVERSE_PW, frameData, &impl_);
    } else {
        const auto encoder = PwEncoder(*this);
        CreateOrEdit(BamRecordTag::REVERSE_PW, encoder.Encode(frameData), &impl_);
    }
    return *this;
}

void BamRecord::ResetCachedPositions() const
{
    alignedEnd_ = Data::UnmappedPosition;
    alignedStart_ = Data::UnmappedPosition;
}

void BamRecord::ResetCachedPositions()
{
    alignedEnd_ = Data::UnmappedPosition;
    alignedStart_ = Data::UnmappedPosition;
}

VirtualRegionType BamRecord::ScrapRegionType() const
{
    const auto tagName = BamRecordTags::LabelFor(BamRecordTag::SCRAP_REGION_TYPE);
    const auto srTag = impl_.TagValue(tagName);
    return VirtualRegionTypeMap::ParseChar[srTag.ToUInt8()];
}

BamRecord& BamRecord::ScrapRegionType(const VirtualRegionType type)
{
    CreateOrEdit(BamRecordTag::SCRAP_REGION_TYPE, static_cast<uint8_t>(type), &impl_);
    return *this;
}

BamRecord& BamRecord::ScrapRegionType(const char type)
{
    CreateOrEdit(BamRecordTag::SCRAP_REGION_TYPE, type, &impl_);
    return *this;
}

ZmwType BamRecord::ScrapZmwType() const
{
    const auto tagName = BamRecordTags::LabelFor(BamRecordTag::SCRAP_ZMW_TYPE);
    const auto szTag = impl_.TagValue(tagName);
    return ZmwTypeMap::ParseChar[szTag.ToUInt8()];
}

BamRecord& BamRecord::ScrapZmwType(const ZmwType type)
{
    CreateOrEdit(BamRecordTag::SCRAP_ZMW_TYPE, static_cast<uint8_t>(type), &impl_);
    return *this;
}

BamRecord& BamRecord::ScrapZmwType(const char type)
{
    CreateOrEdit(BamRecordTag::SCRAP_ZMW_TYPE, type, &impl_);
    return *this;
}

std::string BamRecord::Sequence(const Data::Orientation orientation, bool aligned,
                                bool exciseSoftClips) const
{
    return FetchBases(BamRecordTag::SEQ, orientation, aligned, exciseSoftClips);
}

std::vector<float> BamRecord::SignalToNoise() const
{
    const auto tagName = BamRecordTags::LabelFor(BamRecordTag::SIGNAL_TO_NOISE);
    const auto snTag = impl_.TagValue(tagName);
    return snTag.ToFloatArray();
}

BamRecord& BamRecord::SignalToNoise(const std::vector<float>& snr)
{
    CreateOrEdit(BamRecordTag::SIGNAL_TO_NOISE, snr, &impl_);
    return *this;
}

std::vector<uint32_t> BamRecord::StartFrame(Data::Orientation orientation, bool aligned,
                                            bool exciseSoftClips, PulseBehavior pulseBehavior) const
{
    return FetchUInt32s(BamRecordTag::START_FRAME, orientation, aligned, exciseSoftClips,
                        pulseBehavior);
}

BamRecord& BamRecord::StartFrame(const std::vector<uint32_t>& startFrame)
{
    CreateOrEdit(BamRecordTag::START_FRAME, startFrame, &impl_);
    return *this;
}

QualityValues BamRecord::SubstitutionQV(Data::Orientation orientation, bool aligned,
                                        bool exciseSoftClips) const
{
    return FetchQualities(BamRecordTag::SUBSTITUTION_QV, orientation, aligned, exciseSoftClips);
}

BamRecord& BamRecord::SubstitutionQV(const Data::QualityValues& substitutionQVs)
{
    CreateOrEdit(BamRecordTag::SUBSTITUTION_QV, substitutionQVs.Fastq(), &impl_);
    return *this;
}

std::string BamRecord::SubstitutionTag(Data::Orientation orientation, bool aligned,
                                       bool exciseSoftClips) const
{
    return FetchBases(BamRecordTag::SUBSTITUTION_TAG, orientation, aligned, exciseSoftClips);
}

BamRecord& BamRecord::SubstitutionTag(const std::string& tags)
{
    CreateOrEdit(BamRecordTag::SUBSTITUTION_TAG, tags, &impl_);
    return *this;
}

Data::Read BamRecord::ToRead(std::string model) const
{
    Data::Read result{FullName(),      Sequence(),   Qualities(),
                      SignalToNoise(), QueryStart(), QueryEnd()};
    result.Model = std::move(model);
    result.ReadAccuracy = ReadAccuracy();

    if (HasLocalContextFlags()) {
        result.Flags = LocalContextFlags();
        result.FullLength =
            (result.Flags & Data::ADAPTER_BEFORE) && (result.Flags & Data::ADAPTER_AFTER);
    } else {
        throw std::runtime_error{"[pbbam] BAM record ERROR: '" + FullName() +
                                 "' is missing local context flags (SAM tag 'cx')"};
    }

    if (HasPulseWidth()) {
        result.PulseWidth = PulseWidth();
    }

    if (HasIPD()) {
        result.IPD = IPD();
    }

    if (IsMapped() && AlignedStrand() == Data::Strand::REVERSE) {
        ReverseComplement(result.Seq);
        Reverse(result.Qualities);
    }
    return result;
}

Data::MappedRead BamRecord::ToMappedRead(std::string model, const Data::Position startOffset,
                                         const bool pinStart, const bool pinEnd) const
{
    if (!IsMapped()) {
        throw std::runtime_error{"[pbbam] BAM record ERROR: '" + FullName() +
                                 "' cannot be converted to MappedRead because it is not mapped"};
    }

    Data::MappedRead result{ToRead(std::move(model)), AlignedStrand(), ReferenceStart(),
                            ReferenceEnd(),           CigarData(),     MapQuality()};

    assert(startOffset >= 0);
    result.TemplateStart -= startOffset;
    assert(result.TemplateStart >= 0);
    result.TemplateEnd -= startOffset;
    assert(result.TemplateEnd >= 0);

    result.PinStart = pinStart;
    result.PinEnd = pinEnd;

    return result;
}

RecordType BamRecord::Type() const
{
    try {
        const auto typeName = ReadGroup().ReadType();
        return NameToType(typeName);
    } catch (std::exception&) {

        // read group not found, peek at name to see if we're possibly
        // CCS or TRANSCRIPT
        //
        const auto name = FullName();
        if (name.find("transcript") == 0) {
            return RecordType::TRANSCRIPT;
        } else if (name.find("/ccs") != std::string::npos) {
            return RecordType::CCS;
        } else {
            return RecordType::UNKNOWN;
        }
    }
}

void BamRecord::UpdateName()
{
    std::string newName;
    newName.reserve(100);

    const auto type = Type();
    const auto holeNumber = (HasHoleNumber() ? std::to_string(HoleNumber()) : "?");
    if (type == RecordType::TRANSCRIPT) {
        newName = "transcript/" + holeNumber;
    } else {
        newName += MovieName();
        newName += "/";
        newName += holeNumber;
        newName += "/";

        if (type == RecordType::CCS) {
            newName += "ccs";
            const auto originalName = FullName();
            if (boost::ends_with(originalName, "/fwd")) {
                newName += "/fwd";
            } else if (boost::ends_with(originalName, "/rev")) {
                newName += "/rev";
            }
        }

        else {
            if (HasQueryStart()) {
                newName += std::to_string(QueryStart());
            } else {
                newName += "?";
            }

            newName += '_';

            if (HasQueryEnd()) {
                newName += std::to_string(QueryEnd());
            } else {
                newName += "?";
            }
        }
    }
    impl_.Name(newName);
}

}  // namespace BAM
}  // namespace PacBio
