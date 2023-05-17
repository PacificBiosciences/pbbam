#include "PbbamInternalConfig.h"

#include <pbbam/BamRecordImpl.h>

#include <pbbam/BamTagCodec.h>
#include <pbbam/StringUtilities.h>
#include "BamRecordTags.h"
#include "MemoryUtils.h"

#include <pbcopper/utility/Ssize.h>

#include <htslib/hts_endian.h>

#include <algorithm>
#include <array>
#include <optional>
#include <sstream>
#include <tuple>
#include <utility>

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>

namespace PacBio {
namespace BAM {

namespace {

Data::Cigar FetchRawCigar(const std::uint32_t* const src, const std::uint32_t len)
{
    Data::Cigar result;
    result.reserve(len);
    for (std::uint32_t i = 0; i < len; ++i) {
        const std::uint32_t length = bam_cigar_oplen(src[i]);
        const auto type = static_cast<Data::CigarOperationType>(bam_cigar_op(src[i]));
        result.push_back({type, length});
    }
    return result;
}

}  // namespace

BamRecordImpl::BamRecordImpl() : d_{nullptr}
{
    InitializeData();
    assert(d_);
}

BamRecordImpl::BamRecordImpl(const BamRecordImpl& other)
    : d_{bam_dup1(other.d_.get())}, tagOffsets_{other.tagOffsets_}
{
    assert(d_);
}

BamRecordImpl& BamRecordImpl::operator=(const BamRecordImpl& other)
{
    if (this != &other) {
        if (d_ == nullptr) {
            InitializeData();
        }
        auto* copyOk = bam_copy1(d_.get(), other.d_.get());
        if (!copyOk) {
            throw std::runtime_error{"[pbbam] BAM record ERROR: could not copy data from record '" +
                                     other.Name() + '\''};
        }
        tagOffsets_ = other.tagOffsets_;
    }
    assert(d_);
    return *this;
}

BamRecordImpl::~BamRecordImpl() = default;

bool BamRecordImpl::AddTag(const std::string& tagName, const Tag& value)
{
    return AddTag(tagName, value, TagModifier::NONE);
}

bool BamRecordImpl::AddTag(const BamRecordTag tag, const Tag& value)
{
    return AddTag(BamRecordTags::LabelFor(tag), value, TagModifier::NONE);
}

bool BamRecordImpl::AddTag(const std::string& tagName, const Tag& value,
                           const TagModifier additionalModifier)
{
    if (tagName.size() != 2 || HasTag(tagName)) {
        return false;
    }
    const auto added = AddTagImpl(tagName, value, additionalModifier);
    if (added) {
        UpdateTagMap();
    }
    return added;
}

bool BamRecordImpl::AddTag(const BamRecordTag tag, const Tag& value,
                           const TagModifier additionalModifier)
{
    return AddTag(BamRecordTags::LabelFor(tag), value, additionalModifier);
}

bool BamRecordImpl::AddTagImpl(const std::string& tagName, const Tag& value,
                               const TagModifier additionalModifier)
{
    const auto rawData = BamTagCodec::ToRawData(value, additionalModifier);
    if (rawData.empty()) {
        return false;
    }

    bam_aux_append(d_.get(), tagName.c_str(), BamTagCodec::TagTypeCode(value, additionalModifier),
                   rawData.size(), const_cast<std::uint8_t*>(rawData.data()));
    return true;
}

uint32_t BamRecordImpl::Bin() const { return d_->core.bin; }

BamRecordImpl& BamRecordImpl::Bin(std::uint32_t bin)
{
    d_->core.bin = bin;
    return *this;
}

Data::Cigar BamRecordImpl::CigarData() const
{
    const auto* b = d_.get();
    return FetchRawCigar(bam_get_cigar(b), b->core.n_cigar);
}

BamRecordImpl& BamRecordImpl::CigarData(const Data::Cigar& cigar)
{
    if (HasTag("CG")) {
        RemoveTag("CG");
    }
    SetCigarData(cigar);

    return *this;
}

BamRecordImpl& BamRecordImpl::CigarData(const std::string& cigarString)
{
    return CigarData(Data::Cigar::FromStdString(cigarString));
}

bool BamRecordImpl::EditTag(const std::string& tagName, const Tag& newValue)
{
    return EditTag(tagName, newValue, TagModifier::NONE);
}

bool BamRecordImpl::EditTag(const BamRecordTag tag, const Tag& newValue)
{
    return EditTag(BamRecordTags::LabelFor(tag), newValue, TagModifier::NONE);
}

bool BamRecordImpl::EditTag(const std::string& tagName, const Tag& newValue,
                            const TagModifier additionalModifier)
{
    // try remove old value (with delayed tag map update)
    const bool removed = RemoveTagImpl(tagName);
    if (!removed) {
        return false;
    }

    // if old value removed, add new value
    const bool added = AddTagImpl(tagName, newValue, additionalModifier);
    if (added) {
        UpdateTagMap();
    }
    return added;
}

bool BamRecordImpl::EditTag(const BamRecordTag tag, const Tag& newValue,
                            const TagModifier additionalModifier)
{
    return EditTag(BamRecordTags::LabelFor(tag), newValue, additionalModifier);
}

int BamRecordImpl::EstimatedBytesUsed() const noexcept
{
    // bam1_t
    int result = sizeof(std::unique_ptr<bam1_t, HtslibRecordDeleter>);
    if (d_) {
        bam1_t* b = d_.get();
        result += sizeof(bam1_t);
        result += b->m_data * sizeof(std::uint8_t);  // allocated data block
    }

    // tag offsets cache
    result += sizeof(std::vector<TagOffsetEntry>);
    result += tagOffsets_.capacity() * sizeof(TagOffsetEntry);
    return result;
}

uint32_t BamRecordImpl::Flag() const { return d_->core.flag; }

BamRecordImpl& BamRecordImpl::Flag(std::uint32_t flag)
{
    d_->core.flag = flag;
    return *this;
}

BamRecordImpl BamRecordImpl::FromRawData(const std::shared_ptr<bam1_t>& rawData)
{
    BamRecordImpl result;
    auto* copyOk = bam_copy1(result.d_.get(), rawData.get());
    if (!copyOk) {
        throw std::runtime_error{
            "[pbbam] BAM record ERROR: could not create record, copying from raw BAM contents"};
    }
    return result;
}

bool BamRecordImpl::HasTag(const std::string& tagName) const
{
    if (tagName.size() != 2) {
        return false;
    }
    return TagOffset(tagName) != -1;
}

bool BamRecordImpl::HasTag(const BamRecordTag tag) const
{
    return HasTag(BamRecordTags::LabelFor(tag));
}

void BamRecordImpl::InitializeData()
{
    d_.reset(bam_init1());

    // init unmapped
    Position(Data::UNMAPPED_POSITION);
    MatePosition(Data::UNMAPPED_POSITION);
    ReferenceId(-1);
    MateReferenceId(-1);
    SetMapped(false);
    MapQuality(255);

    // init empty QNAME
    Name("");

    tagOffsets_.reserve(32);
}

int32_t BamRecordImpl::InsertSize() const { return d_->core.isize; }

BamRecordImpl& BamRecordImpl::InsertSize(std::int32_t iSize)
{
    d_->core.isize = iSize;
    return *this;
}

bool BamRecordImpl::IsDuplicate() const { return (d_->core.flag & BamRecordImpl::DUPLICATE) != 0; }

bool BamRecordImpl::IsFailedQC() const { return (d_->core.flag & BamRecordImpl::FAILED_QC) != 0; }

bool BamRecordImpl::IsFirstMate() const { return (d_->core.flag & BamRecordImpl::MATE_1) != 0; }

bool BamRecordImpl::IsMapped() const { return (d_->core.flag & BamRecordImpl::UNMAPPED) == 0; }

bool BamRecordImpl::IsMateMapped() const
{
    return (d_->core.flag & BamRecordImpl::MATE_UNMAPPED) == 0;
}

bool BamRecordImpl::IsMateReverseStrand() const
{
    return (d_->core.flag & BamRecordImpl::MATE_REVERSE_STRAND) != 0;
}

bool BamRecordImpl::IsPaired() const { return (d_->core.flag & BamRecordImpl::PAIRED) != 0; }

bool BamRecordImpl::IsPrimaryAlignment() const
{
    return (d_->core.flag & BamRecordImpl::SECONDARY) == 0;
}

bool BamRecordImpl::IsProperPair() const
{
    return (d_->core.flag & BamRecordImpl::PROPER_PAIR) != 0;
}

bool BamRecordImpl::IsReverseStrand() const
{
    return (d_->core.flag & BamRecordImpl::REVERSE_STRAND) != 0;
}

bool BamRecordImpl::IsSecondMate() const { return (d_->core.flag & BamRecordImpl::MATE_2) != 0; }

bool BamRecordImpl::IsSupplementaryAlignment() const
{
    return (d_->core.flag & BamRecordImpl::SUPPLEMENTARY) != 0;
}

uint8_t BamRecordImpl::MapQuality() const { return d_->core.qual; }

BamRecordImpl& BamRecordImpl::MapQuality(std::uint8_t mapQual)
{
    d_->core.qual = mapQual;
    return *this;
}

Data::Position BamRecordImpl::MatePosition() const { return d_->core.mpos; }

BamRecordImpl& BamRecordImpl::MatePosition(Data::Position pos)
{
    d_->core.mpos = pos;
    return *this;
}

int32_t BamRecordImpl::MateReferenceId() const { return d_->core.mtid; }

BamRecordImpl& BamRecordImpl::MateReferenceId(std::int32_t id)
{
    d_->core.mtid = id;
    return *this;
}

Data::Position BamRecordImpl::Position() const { return d_->core.pos; }

BamRecordImpl& BamRecordImpl::Position(Data::Position pos)
{
    d_->core.pos = pos;
    return *this;
}

int32_t BamRecordImpl::ReferenceId() const { return d_->core.tid; }

BamRecordImpl& BamRecordImpl::ReferenceId(std::int32_t id)
{
    d_->core.tid = id;
    return *this;
}

BamRecordImpl& BamRecordImpl::SetDuplicate(bool ok)
{
    if (ok) {
        d_->core.flag |= BamRecordImpl::DUPLICATE;
    } else {
        d_->core.flag &= ~BamRecordImpl::DUPLICATE;
    }
    return *this;
}

void BamRecordImpl::MaybeReallocData()
{
    // about to grow data contents to l_data size, but m_data is our current max.
    // so we may need to grow. if so, use kroundup to double to next power of 2
    //
    // from sam.h:
    //   decltype(m_data) = uint32_t
    //   decltype(l_data) = int
    if (d_->m_data < static_cast<std::uint32_t>(d_->l_data)) {
        d_->m_data = d_->l_data;
        kroundup32(d_->m_data);
        d_->data = static_cast<std::uint8_t*>(realloc(d_->data, d_->m_data));
    }
}

std::string BamRecordImpl::Name() const { return std::string(bam_get_qname(d_)); }

BamRecordImpl& BamRecordImpl::Name(const std::string& name)
{
    // determine change in memory needed
    // diffNumBytes: pos -> growing, neg -> shrinking
    const auto numChars = Utility::Ssize(name) + 1;  // +1 for NULL-term
    const auto numExtraNulls = 4 - (numChars % 4);
    const auto totalNameSize = numChars + numExtraNulls;

    const auto diffNumBytes = totalNameSize - d_->core.l_qname;
    const auto oldLengthData = d_->l_data;
    d_->l_data += diffNumBytes;
    MaybeReallocData();

    // shift trailing data (cigar, seq, qual, tags) as needed
    const std::uint32_t* oldCigarStart = bam_get_cigar(d_);
    const std::size_t trailingDataLength =
        oldLengthData - (reinterpret_cast<const unsigned char*>(oldCigarStart) -
                         reinterpret_cast<const unsigned char*>(d_->data));
    d_->core.l_qname = totalNameSize;
    d_->core.l_extranul = numExtraNulls;
    std::uint32_t* newCigarStart = bam_get_cigar(d_);
    memmove(newCigarStart, oldCigarStart, trailingDataLength);

    // fill in new name
    memcpy(d_->data, name.c_str(), numChars);
    memset(d_->data + numChars, '\0', numExtraNulls);
    return *this;
}

Data::QualityValues BamRecordImpl::Qualities() const
{
    if (d_->core.l_qseq == 0) {
        return Data::QualityValues();
    }

    std::uint8_t* qualData = bam_get_qual(d_);
    if (qualData[0] == 0xff) {
        return Data::QualityValues();
    }

    const std::size_t numQuals = d_->core.l_qseq;
    Data::QualityValues result;
    result.reserve(numQuals);
    for (std::size_t i = 0; i < numQuals; ++i) {
        result.push_back(Data::QualityValue(qualData[i]));
    }
    return result;
}

bool BamRecordImpl::RemoveTag(const std::string& tagName)
{
    const bool removed = RemoveTagImpl(tagName);
    if (removed) {
        UpdateTagMap();
    }
    return removed;
}

bool BamRecordImpl::RemoveTag(const BamRecordTag tag)
{
    return RemoveTag(BamRecordTags::LabelFor(tag));
}

bool BamRecordImpl::RemoveTagImpl(const std::string& tagName)
{
    if (tagName.size() != 2) {
        return false;
    }
    std::uint8_t* data = bam_aux_get(d_.get(), tagName.c_str());
    if (data == nullptr) {
        return false;
    }
    const bool ok = bam_aux_del(d_.get(), data) == 0;
    return ok;
}

std::string BamRecordImpl::Sequence() const
{
    std::string result(d_->core.l_qseq, '\0');
    constexpr std::array<char, 16> DNA_LOOKUP{
        {'=', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'}};
    const std::uint8_t* seqData = bam_get_seq(d_);
    for (int i = 0; i < d_->core.l_qseq; ++i) {
        result[i] = DNA_LOOKUP[bam_seqi(seqData, i)];
    }
    return result;
}

size_t BamRecordImpl::SequenceLength() const { return d_->core.l_qseq; }

void BamRecordImpl::SetCigarData(const Data::Cigar& cigar)
{
    // determine change in memory needed
    // diffNumBytes: pos -> growing, neg -> shrinking
    const auto numCigarOps = Utility::Ssize(cigar);
    const auto diffNumCigars = numCigarOps - d_->core.n_cigar;
    const auto diffNumBytes = diffNumCigars * static_cast<int>(sizeof(std::uint32_t));
    const auto oldLengthData = d_->l_data;
    d_->l_data += diffNumBytes;
    MaybeReallocData();

    // shift trailing data (seq, qual, tags) as needed
    const std::uint8_t* oldSequenceStart = bam_get_seq(d_);
    const std::size_t trailingDataLength = oldLengthData - (oldSequenceStart - d_->data);
    d_->core.n_cigar = numCigarOps;
    std::uint8_t* newSequenceStart = bam_get_seq(d_);
    memmove(newSequenceStart, oldSequenceStart, trailingDataLength);

    // fill in new CIGAR data
    std::uint32_t* cigarDataStart = bam_get_cigar(d_);
    for (int i = 0; i < numCigarOps; ++i) {
        const Data::CigarOperation& cigarOp = cigar.at(i);
        cigarDataStart[i] = bam_cigar_gen(cigarOp.Length(), static_cast<int>(cigarOp.Type()));
    }
}

BamRecordImpl& BamRecordImpl::SetFailedQC(bool ok)
{
    if (ok) {
        d_->core.flag |= BamRecordImpl::FAILED_QC;
    } else {
        d_->core.flag &= ~BamRecordImpl::FAILED_QC;
    }
    return *this;
}

BamRecordImpl& BamRecordImpl::SetFirstMate(bool ok)
{
    if (ok) {
        d_->core.flag |= BamRecordImpl::MATE_1;
    } else {
        d_->core.flag &= ~BamRecordImpl::MATE_1;
    }
    return *this;
}

BamRecordImpl& BamRecordImpl::SetMapped(bool ok)
{
    if (ok) {
        d_->core.flag &= ~BamRecordImpl::UNMAPPED;
    } else {
        d_->core.flag |= BamRecordImpl::UNMAPPED;
    }
    return *this;
}

BamRecordImpl& BamRecordImpl::SetMateMapped(bool ok)
{
    if (ok) {
        d_->core.flag &= ~BamRecordImpl::MATE_UNMAPPED;
    } else {
        d_->core.flag |= BamRecordImpl::MATE_UNMAPPED;
    }
    return *this;
}

BamRecordImpl& BamRecordImpl::SetMateReverseStrand(bool ok)
{
    if (ok) {
        d_->core.flag |= BamRecordImpl::MATE_REVERSE_STRAND;
    } else {
        d_->core.flag &= ~BamRecordImpl::MATE_REVERSE_STRAND;
    }
    return *this;
}

BamRecordImpl& BamRecordImpl::SetPaired(bool ok)
{
    if (ok) {
        d_->core.flag |= BamRecordImpl::PAIRED;
    } else {
        d_->core.flag &= ~BamRecordImpl::PAIRED;
    }
    return *this;
}

BamRecordImpl& BamRecordImpl::SetPrimaryAlignment(bool ok)
{
    if (ok) {
        d_->core.flag &= ~BamRecordImpl::SECONDARY;
    } else {
        d_->core.flag |= BamRecordImpl::SECONDARY;
    }
    return *this;
}

BamRecordImpl& BamRecordImpl::SetProperPair(bool ok)
{
    if (ok) {
        d_->core.flag |= BamRecordImpl::PROPER_PAIR;
    } else {
        d_->core.flag &= ~BamRecordImpl::PROPER_PAIR;
    }
    return *this;
}

BamRecordImpl& BamRecordImpl::SetReverseStrand(bool ok)
{
    if (ok) {
        d_->core.flag |= BamRecordImpl::REVERSE_STRAND;
    } else {
        d_->core.flag &= ~BamRecordImpl::REVERSE_STRAND;
    }
    return *this;
}

BamRecordImpl& BamRecordImpl::SetSecondMate(bool ok)
{
    if (ok) {
        d_->core.flag |= BamRecordImpl::MATE_2;
    } else {
        d_->core.flag &= ~BamRecordImpl::MATE_2;
    }
    return *this;
}

BamRecordImpl& BamRecordImpl::SetSequenceAndQualities(const std::string& sequence,
                                                      const std::string& qualities)
{
    if (!qualities.empty() && (sequence.size() != qualities.size())) {
        std::ostringstream s;
        s << "[pbbam] BAM record ERROR: if qualities are provided, the length must match the "
             "sequence length:\n"
          << "  sequence length: " << sequence.size() << '\n'
          << "  qualities length: " << qualities.size();
        throw std::runtime_error{s.str()};
    }
    return SetSequenceAndQualitiesInternal(sequence.c_str(), sequence.size(), qualities.c_str(),
                                           false);
}

BamRecordImpl& BamRecordImpl::SetSequenceAndQualities(const char* sequence,
                                                      const std::size_t sequenceLength,
                                                      const char* qualities)
{
    return SetSequenceAndQualitiesInternal(sequence, sequenceLength, qualities, false);
}

BamRecordImpl& BamRecordImpl::SetPreencodedSequenceAndQualities(const char* encodedSequence,
                                                                const std::size_t rawSequenceLength,
                                                                const char* qualities)
{
    return SetSequenceAndQualitiesInternal(encodedSequence, rawSequenceLength, qualities, true);
}

BamRecordImpl& BamRecordImpl::SetSequenceAndQualitiesInternal(const char* sequence,
                                                              const std::size_t sequenceLength,
                                                              const char* qualities,
                                                              bool isPreencoded)
{
    // determine change in memory needed
    // diffNumBytes: pos -> growing, neg -> shrinking
    const auto encodedSequenceLength = static_cast<int>((sequenceLength + 1) / 2);
    const int oldSeqAndQualLength =
        ((d_->core.l_qseq + 1) / 2) + d_->core.l_qseq;                       // encoded seq + qual
    const int newSeqAndQualLength = encodedSequenceLength + sequenceLength;  // encoded seq + qual
    const int diffNumBytes = newSeqAndQualLength - oldSeqAndQualLength;
    const int oldLengthData = d_->l_data;
    d_->l_data += diffNumBytes;
    MaybeReallocData();

    // shift trailing data (tags) as needed
    const unsigned char* oldTagStart = bam_get_aux(d_);
    const std::size_t trailingDataLength =
        oldLengthData - (oldTagStart - reinterpret_cast<const unsigned char*>(d_->data));
    d_->core.l_qseq = sequenceLength;
    std::uint8_t* newTagStart = bam_get_aux(d_);
    memmove(newTagStart, oldTagStart, trailingDataLength);

    // fill in new sequence
    std::uint8_t* pEncodedSequence = bam_get_seq(d_);
    if (isPreencoded) {
        memcpy(pEncodedSequence, sequence, encodedSequenceLength);
    } else {
        memset(pEncodedSequence, 0, encodedSequenceLength);
        for (std::size_t i = 0; i < sequenceLength; ++i) {
            pEncodedSequence[i >> 1] |= seq_nt16_table[static_cast<int>(sequence[i])]
                                        << ((~i & 1) << 2);
        }
    }

    // fill in quality values
    std::uint8_t* encodedQualities = bam_get_qual(d_);
    if ((qualities == nullptr) || (strlen(qualities) == 0)) {
        memset(encodedQualities, 0xff, sequenceLength);
    } else {
        for (std::size_t i = 0; i < sequenceLength; ++i) {
            encodedQualities[i] = qualities[i] - 33;  // FASTQ ASCII -> int conversion
        }
    }
    return *this;
}

BamRecordImpl& BamRecordImpl::SetSupplementaryAlignment(bool ok)
{
    if (ok) {
        d_->core.flag |= BamRecordImpl::SUPPLEMENTARY;
    } else {
        d_->core.flag &= ~BamRecordImpl::SUPPLEMENTARY;
    }
    return *this;
}

std::optional<int> BamRecordImpl::TagLength(const std::string& tagName) const
{
    const int offset = TagOffset(tagName);

    // not present
    if (offset == -1) {
        return std::nullopt;
    }

    // fetch tag data section
    bam1_t* b = d_.get();
    assert(bam_get_aux(b));
    std::uint8_t* tagData = bam_get_aux(b) + offset;
    if (offset >= b->l_data) {
        return std::nullopt;
    }

    // determine tag length
    const auto tagType = static_cast<char>(*tagData++);
    switch (tagType) {

        // string tag
        case 'H':
            [[fallthrough]];
        case 'Z': {
            return strlen(reinterpret_cast<const char*>(&tagData[0]));
        }

        // array tag
        case 'B': {
            ++tagData;  // skip array's value type code
            std::uint32_t numElements = 0;
            memcpy(&numElements, &tagData[0], sizeof(std::uint32_t));
            return numElements;
        }

        // scalar tag
        default:
            return std::nullopt;
    }
}

std::optional<int> BamRecordImpl::TagLength(const BamRecordTag tag) const
{
    return TagLength(BamRecordTags::LabelFor(tag));
}

int BamRecordImpl::TagOffset(const std::string& tagName) const
{
    if (tagName.size() != 2) {
        throw std::runtime_error{"[pbbam] BAM record ERROR: tag name (" + tagName +
                                 ") must have 2 characters only"};
    }

    if (tagOffsets_.empty()) {
        UpdateTagMap();
    }

    const std::uint16_t tagCode =
        (static_cast<std::uint8_t>(tagName.at(0)) << 8) | static_cast<std::uint8_t>(tagName.at(1));
    for (const auto& tag : tagOffsets_) {
        if (tag.Code == tagCode) {
            return tag.Offset;
        }
    }
    return -1;  // not found
}

BamRecordImpl& BamRecordImpl::Tags(const TagCollection& tags)
{
    // convert tags to binary
    const std::vector<std::uint8_t> tagData = BamTagCodec::Encode(tags);
    const int numBytes = Utility::Ssize(tagData);
    const std::uint8_t* data = tagData.data();

    // determine change in memory needed
    std::uint8_t* tagStart = bam_get_aux(d_);
    const int oldNumBytes = d_->l_data - (tagStart - d_->data);
    const int diffNumBytes = numBytes - oldNumBytes;
    d_->l_data += diffNumBytes;
    MaybeReallocData();
    tagStart = bam_get_aux(d_);

    // fill in new tag data
    if (numBytes) {
        std::memcpy(static_cast<void*>(tagStart), data, numBytes);
    }

    // update tag info
    UpdateTagMap();
    return *this;
}

TagCollection BamRecordImpl::Tags() const
{
    const std::uint8_t* tagDataStart = bam_get_aux(d_);
    const std::size_t numBytes = d_->l_data - (tagDataStart - d_->data);
    return BamTagCodec::Decode(std::vector<std::uint8_t>(tagDataStart, tagDataStart + numBytes));
}

Tag BamRecordImpl::TagValue(const std::string& tagName) const
{
    if (tagName.size() != 2) {
        return {};
    }

    const int offset = TagOffset(tagName);
    if (offset == -1) {
        return {};
    }

    bam1_t* b = d_.get();
    assert(bam_get_aux(b));
    std::uint8_t* tagData = bam_get_aux(b) + offset;
    if (offset >= b->l_data) {
        return {};
    }

    // skip tag name
    return BamTagCodec::FromRawData(tagData);
}

Tag BamRecordImpl::TagValue(const BamRecordTag tag) const
{
    return TagValue(BamRecordTags::LabelFor(tag));
}

void BamRecordImpl::UpdateTagMap() const
{
    tagOffsets_.clear();

    const std::uint8_t* tagStart = bam_get_aux(d_);
    if (tagStart == nullptr) {
        return;
    }
    const std::ptrdiff_t numBytes = d_->l_data - (tagStart - d_->data);

    // NOTE: using a 16-bit 'code' for tag name here instead of string, to avoid
    // a lot of string constructions & comparisons. All valid tags will be 2 chars
    // anyway, so this should be a nice lookup mechanism.
    //
    std::uint16_t tagNameCode;
    int i = 0;
    while (i < numBytes) {

        // store (tag name code -> start offset into tag data)
        tagNameCode = static_cast<char>(tagStart[i]) << 8 | static_cast<char>(tagStart[i + 1]);
        i += 2;
        tagOffsets_.push_back(TagOffsetEntry{tagNameCode, i});

        // skip tag contents
        const auto tagType = static_cast<char>(tagStart[i++]);
        switch (tagType) {
            case 'A':
            case 'a':
            case 'c':
            case 'C': {
                i += 1;
                break;
            }
            case 's':
            case 'S': {
                i += 2;
                break;
            }
            case 'i':
            case 'I':
            case 'f': {
                i += 4;
                break;
            }

            case 'Z':
            case 'H': {
                // null-terminated string
                i += strlen(reinterpret_cast<const char*>(&tagStart[i])) + 1;
                break;
            }

            case 'B': {
                const char subTagType = tagStart[i++];
                std::size_t elementSize = 0;
                switch (subTagType) {
                    case 'c':
                    case 'C':
                        elementSize = 1;
                        break;
                    case 's':
                    case 'S':
                        elementSize = 2;
                        break;
                    case 'i':
                    case 'I':
                    case 'f':
                        elementSize = 4;
                        break;

                    // unknown subTagType
                    default:
                        throw std::runtime_error{
                            "[pbbam] BAM record ERROR: unsupported array-tag-type encountered: " +
                            std::string{1, subTagType}};
                }

                std::uint32_t numElements = 0;
                memcpy(&numElements, &tagStart[i], sizeof(std::uint32_t));
                i += (4 + (elementSize * numElements));
                break;
            }

            // unknown tagType
            default:
                throw std::runtime_error{
                    "[pbbam] BAM record ERROR: unsupported tag-type encountered: " +
                    std::string{1, tagType}};
        }
    }
}

}  // namespace BAM
}  // namespace PacBio
