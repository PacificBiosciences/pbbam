#include "PbbamInternalConfig.h"

#include <pbbam/BamRecordImpl.h>

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <array>
#include <sstream>
#include <tuple>
#include <utility>

#include <htslib/hts_endian.h>

#include <pbbam/BamTagCodec.h>
#include <pbbam/StringUtilities.h>

#include "BamRecordTags.h"
#include "MemoryUtils.h"

namespace PacBio {
namespace BAM {

namespace {

static const bool has_native_long_cigar_support = DoesHtslibSupportLongCigar();

Data::Cigar FetchRawCigar(const uint32_t* const src, const uint32_t len)
{
    Data::Cigar result;
    result.reserve(len);
    for (uint32_t i = 0; i < len; ++i) {
        const uint32_t length = bam_cigar_oplen(src[i]);
        const auto type = static_cast<Data::CigarOperationType>(bam_cigar_op(src[i]));
        result.push_back({type, length});
    }
    return result;
}

bool HasLongCigar(const bam1_t* const b)
{
    auto* c = &b->core;

    // if empty CIGAR or unmapped
    if (c->n_cigar == 0 || c->tid < 0 || c->pos < 0) return false;

    // if existing CIGAR doesn't look like a 'fake CIGAR'
    const auto firstCigarOp = *(bam_get_cigar(b));
    if (bam_cigar_op(firstCigarOp) != static_cast<uint32_t>(Data::CigarOperationType::SOFT_CLIP) ||
        static_cast<int32_t>(bam_cigar_oplen(firstCigarOp)) != c->l_qseq) {
        return false;
    }

    // if CG tag missing, not expected type
    const uint8_t* const CG = bam_aux_get(b, "CG");
    if (CG == nullptr) return false;
    if (CG[0] != 'B' || CG[1] != 'I') return false;

    // if CG tag data is empty
    uint32_t numElements = 0;
    memcpy(&numElements, &CG[2], sizeof(uint32_t));
    if (numElements == 0) return false;

    // we've found long CIGAR data in the CG tag
    return true;
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
        if (d_ == nullptr) InitializeData();
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
    if (tagName.size() != 2 || HasTag(tagName)) return false;
    const auto added = AddTagImpl(tagName, value, additionalModifier);
    if (added) UpdateTagMap();
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
    if (rawData.empty()) return false;

    bam_aux_append(d_.get(), tagName.c_str(), BamTagCodec::TagTypeCode(value, additionalModifier),
                   rawData.size(), const_cast<uint8_t*>(rawData.data()));
    return true;
}

uint32_t BamRecordImpl::Bin() const { return d_->core.bin; }

BamRecordImpl& BamRecordImpl::Bin(uint32_t bin)
{
    d_->core.bin = bin;
    return *this;
}

Data::Cigar BamRecordImpl::CigarData() const
{
    const auto* b = d_.get();
    if (!has_native_long_cigar_support && HasLongCigar(b)) {
        // fetch long CIGAR from tag
        const auto cigarTag = TagValue("CG");
        const auto cigarTagValue = cigarTag.ToUInt32Array();
        return FetchRawCigar(cigarTagValue.data(), cigarTagValue.size());
    } else {
        // fetch CIGAR from the standard location
        return FetchRawCigar(bam_get_cigar(b), b->core.n_cigar);
    }
}

BamRecordImpl& BamRecordImpl::CigarData(const Data::Cigar& cigar)
{
    // if long CIGAR, using htslib version < 1.7, set it "manually"
    if (!has_native_long_cigar_support && cigar.size() >= 65536) {
        // Add the 'fake' CIGAR in normal place.
        Data::Cigar fake;
        fake.emplace_back(Data::CigarOperationType::SOFT_CLIP, SequenceLength());
        const uint32_t alignedLength =
            static_cast<uint32_t>(bam_cigar2rlen(d_->core.n_cigar, bam_get_cigar(d_.get())));
        fake.emplace_back(Data::CigarOperationType::REFERENCE_SKIP, alignedLength);
        SetCigarData(fake);

        // Add raw CIGAR data to CG tag.
        std::vector<uint32_t> cigarData(cigar.size());
        cigarData.reserve(cigar.size());
        for (size_t i = 0; i < cigar.size(); ++i) {
            const CigarOperation& op = cigar.at(i);
            cigarData[i] = bam_cigar_gen(op.Length(), static_cast<int>(op.Type()));
        }
        if (HasTag("CG"))
            EditTag("CG", Tag{cigarData});
        else
            AddTag("CG", Tag{cigarData});
    }

    // otherwise (v1.7+ or short CIGAR), use standard APIs
    else {
        if (HasTag("CG")) RemoveTag("CG");
        SetCigarData(cigar);
    }

    return *this;
}

BamRecordImpl& BamRecordImpl::CigarData(const std::string& cigarString)
{
    return CigarData(Cigar::FromStdString(cigarString));
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
    if (!removed) return false;

    // if old value removed, add new value
    const bool added = AddTagImpl(tagName, newValue, additionalModifier);
    if (added) UpdateTagMap();
    return added;
}

bool BamRecordImpl::EditTag(const BamRecordTag tag, const Tag& newValue,
                            const TagModifier additionalModifier)
{
    return EditTag(BamRecordTags::LabelFor(tag), newValue, additionalModifier);
}

uint32_t BamRecordImpl::Flag() const { return d_->core.flag; }

BamRecordImpl& BamRecordImpl::Flag(uint32_t flag)
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
    if (tagName.size() != 2) return false;
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
    Position(Data::UnmappedPosition);
    MatePosition(Data::UnmappedPosition);
    ReferenceId(-1);
    MateReferenceId(-1);
    SetMapped(false);
    MapQuality(255);

    // init empty QNAME
    Name("");
}

int32_t BamRecordImpl::InsertSize() const { return d_->core.isize; }

BamRecordImpl& BamRecordImpl::InsertSize(int32_t iSize)
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

BamRecordImpl& BamRecordImpl::MapQuality(uint8_t mapQual)
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

BamRecordImpl& BamRecordImpl::MateReferenceId(int32_t id)
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

BamRecordImpl& BamRecordImpl::ReferenceId(int32_t id)
{
    d_->core.tid = id;
    return *this;
}

BamRecordImpl& BamRecordImpl::SetDuplicate(bool ok)
{
    if (ok)
        d_->core.flag |= BamRecordImpl::DUPLICATE;
    else
        d_->core.flag &= ~BamRecordImpl::DUPLICATE;
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
    if (d_->m_data < static_cast<uint32_t>(d_->l_data)) {
        d_->m_data = d_->l_data;
        kroundup32(d_->m_data);
        d_->data = static_cast<uint8_t*>(realloc(d_->data, d_->m_data));
    }
}

std::string BamRecordImpl::Name() const { return std::string(bam_get_qname(d_)); }

BamRecordImpl& BamRecordImpl::Name(const std::string& name)
{
    // determine change in memory needed
    // diffNumBytes: pos -> growing, neg -> shrinking
    const size_t numChars = name.size() + 1;  // +1 for NULL-term
    const size_t numExtraNulls = 4 - (numChars % 4);
    const size_t totalNameSize = numChars + numExtraNulls;

    const int diffNumBytes = totalNameSize - d_->core.l_qname;
    const int oldLengthData = d_->l_data;
    d_->l_data += diffNumBytes;
    MaybeReallocData();

    // shift trailing data (cigar, seq, qual, tags) as needed
    const uint32_t* oldCigarStart = bam_get_cigar(d_);
    const size_t trailingDataLength =
        oldLengthData - (reinterpret_cast<const unsigned char*>(oldCigarStart) -
                         reinterpret_cast<const unsigned char*>(d_->data));
    d_->core.l_qname = totalNameSize;
    d_->core.l_extranul = numExtraNulls;
    uint32_t* newCigarStart = bam_get_cigar(d_);
    memmove(newCigarStart, oldCigarStart, trailingDataLength);

    // fill in new name
    memcpy(d_->data, name.c_str(), numChars);
    memset(d_->data + numChars, '\0', numExtraNulls);
    return *this;
}

Data::QualityValues BamRecordImpl::Qualities() const
{
    if (d_->core.l_qseq == 0) return Data::QualityValues();

    uint8_t* qualData = bam_get_qual(d_);
    if (qualData[0] == 0xff) return Data::QualityValues();

    const size_t numQuals = d_->core.l_qseq;
    Data::QualityValues result;
    result.reserve(numQuals);
    for (size_t i = 0; i < numQuals; ++i)
        result.push_back(Data::QualityValue(qualData[i]));
    return result;
}

bool BamRecordImpl::RemoveTag(const std::string& tagName)
{
    const bool removed = RemoveTagImpl(tagName);
    if (removed) UpdateTagMap();
    return removed;
}

bool BamRecordImpl::RemoveTag(const BamRecordTag tag)
{
    return RemoveTag(BamRecordTags::LabelFor(tag));
}

bool BamRecordImpl::RemoveTagImpl(const std::string& tagName)
{
    if (tagName.size() != 2) return false;
    uint8_t* data = bam_aux_get(d_.get(), tagName.c_str());
    if (data == nullptr) return false;
    const bool ok = bam_aux_del(d_.get(), data) == 0;
    return ok;
}

std::string BamRecordImpl::Sequence() const
{
    std::string result(d_->core.l_qseq, '\0');
    static const constexpr std::array<char, 16> DnaLookup{
        {'=', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'}};
    const uint8_t* seqData = bam_get_seq(d_);
    for (int i = 0; i < d_->core.l_qseq; ++i)
        result[i] = DnaLookup[bam_seqi(seqData, i)];
    return result;
}

size_t BamRecordImpl::SequenceLength() const { return d_->core.l_qseq; }

void BamRecordImpl::SetCigarData(const Data::Cigar& cigar)
{
    // determine change in memory needed
    // diffNumBytes: pos -> growing, neg -> shrinking
    const size_t numCigarOps = cigar.size();
    const int diffNumCigars = numCigarOps - d_->core.n_cigar;
    const int diffNumBytes = diffNumCigars * sizeof(uint32_t);
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
        const Data::CigarOperation& cigarOp = cigar.at(i);
        cigarDataStart[i] = bam_cigar_gen(cigarOp.Length(), static_cast<int>(cigarOp.Type()));
    }
}

BamRecordImpl& BamRecordImpl::SetFailedQC(bool ok)
{
    if (ok)
        d_->core.flag |= BamRecordImpl::FAILED_QC;
    else
        d_->core.flag &= ~BamRecordImpl::FAILED_QC;
    return *this;
}

BamRecordImpl& BamRecordImpl::SetFirstMate(bool ok)
{
    if (ok)
        d_->core.flag |= BamRecordImpl::MATE_1;
    else
        d_->core.flag &= ~BamRecordImpl::MATE_1;
    return *this;
}

BamRecordImpl& BamRecordImpl::SetMapped(bool ok)
{
    if (ok)
        d_->core.flag &= ~BamRecordImpl::UNMAPPED;
    else
        d_->core.flag |= BamRecordImpl::UNMAPPED;
    return *this;
}

BamRecordImpl& BamRecordImpl::SetMateMapped(bool ok)
{
    if (ok)
        d_->core.flag &= ~BamRecordImpl::MATE_UNMAPPED;
    else
        d_->core.flag |= BamRecordImpl::MATE_UNMAPPED;
    return *this;
}

BamRecordImpl& BamRecordImpl::SetMateReverseStrand(bool ok)
{
    if (ok)
        d_->core.flag |= BamRecordImpl::MATE_REVERSE_STRAND;
    else
        d_->core.flag &= ~BamRecordImpl::MATE_REVERSE_STRAND;
    return *this;
}

BamRecordImpl& BamRecordImpl::SetPaired(bool ok)
{
    if (ok)
        d_->core.flag |= BamRecordImpl::PAIRED;
    else
        d_->core.flag &= ~BamRecordImpl::PAIRED;
    return *this;
}

BamRecordImpl& BamRecordImpl::SetPrimaryAlignment(bool ok)
{
    if (ok)
        d_->core.flag &= ~BamRecordImpl::SECONDARY;
    else
        d_->core.flag |= BamRecordImpl::SECONDARY;
    return *this;
}

BamRecordImpl& BamRecordImpl::SetProperPair(bool ok)
{
    if (ok)
        d_->core.flag |= BamRecordImpl::PROPER_PAIR;
    else
        d_->core.flag &= ~BamRecordImpl::PROPER_PAIR;
    return *this;
}

BamRecordImpl& BamRecordImpl::SetReverseStrand(bool ok)
{
    if (ok)
        d_->core.flag |= BamRecordImpl::REVERSE_STRAND;
    else
        d_->core.flag &= ~BamRecordImpl::REVERSE_STRAND;
    return *this;
}

BamRecordImpl& BamRecordImpl::SetSecondMate(bool ok)
{
    if (ok)
        d_->core.flag |= BamRecordImpl::MATE_2;
    else
        d_->core.flag &= ~BamRecordImpl::MATE_2;
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
                                                      const size_t sequenceLength,
                                                      const char* qualities)
{
    return SetSequenceAndQualitiesInternal(sequence, sequenceLength, qualities, false);
}

BamRecordImpl& BamRecordImpl::SetPreencodedSequenceAndQualities(const char* encodedSequence,
                                                                const size_t rawSequenceLength,
                                                                const char* qualities)
{
    return SetSequenceAndQualitiesInternal(encodedSequence, rawSequenceLength, qualities, true);
}

BamRecordImpl& BamRecordImpl::SetSequenceAndQualitiesInternal(const char* sequence,
                                                              const size_t sequenceLength,
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
    const size_t trailingDataLength =
        oldLengthData - (oldTagStart - reinterpret_cast<const unsigned char*>(d_->data));
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
            pEncodedSequence[i >> 1] |= seq_nt16_table[static_cast<int>(sequence[i])]
                                        << ((~i & 1) << 2);
    }

    // fill in quality values
    uint8_t* encodedQualities = bam_get_qual(d_);
    if ((qualities == nullptr) || (strlen(qualities) == 0))
        memset(encodedQualities, 0xff, sequenceLength);
    else {
        for (size_t i = 0; i < sequenceLength; ++i)
            encodedQualities[i] = qualities[i] - 33;  // FASTQ ASCII -> int conversion
    }
    return *this;
}

BamRecordImpl& BamRecordImpl::SetSupplementaryAlignment(bool ok)
{
    if (ok)
        d_->core.flag |= BamRecordImpl::SUPPLEMENTARY;
    else
        d_->core.flag &= ~BamRecordImpl::SUPPLEMENTARY;
    return *this;
}

int BamRecordImpl::TagOffset(const std::string& tagName) const
{
    if (tagName.size() != 2)
        throw std::runtime_error{"[pbbam] BAM record ERROR: tag name (" + tagName +
                                 ") must have 2 characters only"};

    if (tagOffsets_.empty()) UpdateTagMap();

    const uint16_t tagCode =
        (static_cast<uint8_t>(tagName.at(0)) << 8) | static_cast<uint8_t>(tagName.at(1));
    const auto found = tagOffsets_.find(tagCode);
    return (found != tagOffsets_.cend() ? found->second : -1);
}

BamRecordImpl& BamRecordImpl::Tags(const TagCollection& tags)
{
    // convert tags to binary
    const std::vector<uint8_t> tagData = BamTagCodec::Encode(tags);
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
    if (numBytes) {
        std::memcpy(static_cast<void*>(tagStart), data, numBytes);
    }

    // update tag info
    UpdateTagMap();
    return *this;
}

TagCollection BamRecordImpl::Tags() const
{
    const uint8_t* tagDataStart = bam_get_aux(d_);
    const size_t numBytes = d_->l_data - (tagDataStart - d_->data);
    return BamTagCodec::Decode(std::vector<uint8_t>(tagDataStart, tagDataStart + numBytes));
}

Tag BamRecordImpl::TagValue(const std::string& tagName) const
{
    if (tagName.size() != 2) return {};

    const int offset = TagOffset(tagName);
    if (offset == -1) return {};

    bam1_t* b = d_.get();
    assert(bam_get_aux(b));
    uint8_t* tagData = bam_get_aux(b) + offset;
    if (offset >= b->l_data) return {};

    // skip tag name
    return BamTagCodec::FromRawData(tagData);
}

Tag BamRecordImpl::TagValue(const BamRecordTag tag) const
{
    return TagValue(BamRecordTags::LabelFor(tag));
}

void BamRecordImpl::UpdateTagMap() const
{
    // clear out offsets, leave map structure basically intact
    for (auto& tag : tagOffsets_)
        tag.second = -1;

    const uint8_t* tagStart = bam_get_aux(d_);
    if (tagStart == nullptr) return;
    const ptrdiff_t numBytes = d_->l_data - (tagStart - d_->data);

    // NOTE: using a 16-bit 'code' for tag name here instead of string, to avoid
    // a lot of string constructions & comparisons. All valid tags will be 2 chars
    // anyway, so this should be a nice lookup mechanism.
    //
    uint16_t tagNameCode;
    int64_t i = 0;
    while (i < numBytes) {

        // store (tag name code -> start offset into tag data)
        tagNameCode = static_cast<char>(tagStart[i]) << 8 | static_cast<char>(tagStart[i + 1]);
        i += 2;
        tagOffsets_[tagNameCode] = i;

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
                size_t elementSize = 0;
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

                uint32_t numElements = 0;
                memcpy(&numElements, &tagStart[i], sizeof(uint32_t));
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
