// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/BamRecordBuilder.h"

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>

#include <htslib/sam.h>

#include "MemoryUtils.h"
#include "pbbam/BamTagCodec.h"

namespace PacBio {
namespace BAM {

BamRecordBuilder::BamRecordBuilder()
{
    // ensure proper clean slate
    Reset();

    // initialize with some space for data
    name_.reserve(256);
    sequence_.reserve(2096);
    qualities_.reserve(2096);
    cigar_.reserve(256);
}

BamRecordBuilder::BamRecordBuilder(BamHeader header) : header_{std::move(header)}
{
    // ensure proper clean slate
    Reset();

    // initialize with some space for data
    name_.reserve(256);
    sequence_.reserve(2096);
    qualities_.reserve(2096);
    cigar_.reserve(256);
}

BamRecordBuilder::BamRecordBuilder(const BamRecord& prototype) : header_{prototype.Header()}
{
    Reset(prototype);
}

BamRecord BamRecordBuilder::Build() const
{
    BamRecord result{header_};
    BuildInPlace(result);
    return result;
}

bool BamRecordBuilder::BuildInPlace(BamRecord& record) const
{
    // initialize with basic 'core data'
    auto recordRawData = internal::BamRecordMemory::GetRawData(record);
    if (!recordRawData || !recordRawData->data)
        throw std::runtime_error{"BamRecord memory in invalid state"};
    recordRawData->core = core_;

    // setup variable length data
    const auto encodedTags = BamTagCodec::Encode(tags_);

    const size_t nameLength = name_.size() + 1;
    const size_t numCigarOps = cigar_.size();
    const size_t cigarLength = numCigarOps * sizeof(uint32_t);
    const size_t seqLength = sequence_.size();
    const size_t qualLength = seqLength;
    const size_t tagLength = encodedTags.size();
    const size_t dataLength = nameLength + cigarLength + seqLength + qualLength + tagLength;

    // realloc if necessary
    uint8_t* varLengthDataBlock = recordRawData->data;
    if (!varLengthDataBlock) throw std::runtime_error{"BamRecord memory in invalid state"};

    size_t allocatedDataLength = recordRawData->m_data;
    if (allocatedDataLength < dataLength) {
        allocatedDataLength = dataLength;
        kroundup32(allocatedDataLength);
        varLengthDataBlock =
            static_cast<uint8_t*>(realloc(varLengthDataBlock, allocatedDataLength));
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
        std::vector<uint32_t> encodedCigar(numCigarOps);
        for (size_t i = 0; i < numCigarOps; ++i) {
            const auto& op = cigar_.at(i);
            encodedCigar[i] = op.Length() << BAM_CIGAR_SHIFT;
            const auto type = static_cast<uint8_t>(op.Type());
            if (type >= 8)
                throw std::runtime_error{"invalid CIGAR op type: " + std::to_string(type)};
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
            s[i >> 1] |= (seq_nt16_table[static_cast<int>(sequence_.at(i))] << ((~i & 1) << 2));
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
        if (encodedTags.empty()) throw std::runtime_error{"expected tags but none are encoded"};
        memcpy(&varLengthDataBlock[index], &encodedTags[0], tagLength);
        index += tagLength;
    }

    // sanity check
    if (index != dataLength) {
        throw std::runtime_error{"BAM encoding error: expected to write " +
                                 std::to_string(dataLength) + " bytes but wrote " +
                                 std::to_string(index) + " bytes instead"};
    }
    return true;
}

BamRecordBuilder& BamRecordBuilder::Cigar(PacBio::BAM::Cigar cigar)
{
    core_.n_cigar = cigar.size();
    cigar_ = std::move(cigar);
    return *this;
}

BamRecordBuilder& BamRecordBuilder::Name(std::string name)
{
    core_.l_qname = name.size() + 1;  // (NULL-term)
    name_ = std::move(name);
    return *this;
}

void BamRecordBuilder::Reset()
{
    // zeroize fixed-length data
    memset(&core_, 0, sizeof(bam1_core_t));
    core_.l_qname = 1;  // always has a NULL-term

    // reset variable-length data
    name_.clear();
    sequence_.clear();
    qualities_.clear();
    cigar_.clear();
    tags_.clear();
}

void BamRecordBuilder::Reset(BamRecord prototype)
{
    // ensure clean slate
    Reset();
    header_ = prototype.Header();

    // reset variable-length data
    const BamRecordImpl& impl = internal::BamRecordMemory::GetImpl(prototype);
    name_ = impl.Name();
    sequence_ = impl.Sequence();
    qualities_ = impl.Qualities().Fastq();
    cigar_ = impl.CigarData();
    tags_ = impl.Tags();

    // reset core data
    const auto rawData =
        internal::BamRecordMemory::GetRawData(prototype);  //  prototype.impl_.RawData().get();
    if (!rawData) throw std::runtime_error{"BamRecord memory in invalid state"};
    core_ = std::move(rawData->core);
}

BamRecordBuilder& BamRecordBuilder::Sequence(std::string sequence)
{
    core_.l_qseq = sequence.size();
    sequence_ = std::move(sequence);
    return *this;
}

BamRecordBuilder& BamRecordBuilder::SetDuplicate(bool ok)
{
    if (ok)
        core_.flag |= BamRecordImpl::DUPLICATE;
    else
        core_.flag &= ~BamRecordImpl::DUPLICATE;
    return *this;
}

BamRecordBuilder& BamRecordBuilder::SetFailedQC(bool ok)
{
    if (ok)
        core_.flag |= BamRecordImpl::FAILED_QC;
    else
        core_.flag &= ~BamRecordImpl::FAILED_QC;
    return *this;
}

BamRecordBuilder& BamRecordBuilder::SetFirstMate(bool ok)
{
    if (ok)
        core_.flag |= BamRecordImpl::MATE_1;
    else
        core_.flag &= ~BamRecordImpl::MATE_1;
    return *this;
}

BamRecordBuilder& BamRecordBuilder::SetMapped(bool ok)
{
    if (ok)
        core_.flag &= ~BamRecordImpl::UNMAPPED;
    else
        core_.flag |= BamRecordImpl::UNMAPPED;
    return *this;
}

BamRecordBuilder& BamRecordBuilder::SetMateMapped(bool ok)
{
    if (ok)
        core_.flag &= ~BamRecordImpl::MATE_UNMAPPED;
    else
        core_.flag |= BamRecordImpl::MATE_UNMAPPED;
    return *this;
}

BamRecordBuilder& BamRecordBuilder::SetMateReverseStrand(bool ok)
{
    if (ok)
        core_.flag |= BamRecordImpl::MATE_REVERSE_STRAND;
    else
        core_.flag &= ~BamRecordImpl::MATE_REVERSE_STRAND;
    return *this;
}

BamRecordBuilder& BamRecordBuilder::SetPaired(bool ok)
{
    if (ok)
        core_.flag |= BamRecordImpl::PAIRED;
    else
        core_.flag &= ~BamRecordImpl::PAIRED;
    return *this;
}

BamRecordBuilder& BamRecordBuilder::SetPrimaryAlignment(bool ok)
{
    if (ok)
        core_.flag &= ~BamRecordImpl::SECONDARY;
    else
        core_.flag |= BamRecordImpl::SECONDARY;
    return *this;
}

BamRecordBuilder& BamRecordBuilder::SetProperPair(bool ok)
{
    if (ok)
        core_.flag |= BamRecordImpl::PROPER_PAIR;
    else
        core_.flag &= ~BamRecordImpl::PROPER_PAIR;
    return *this;
}

BamRecordBuilder& BamRecordBuilder::SetReverseStrand(bool ok)
{
    if (ok)
        core_.flag |= BamRecordImpl::REVERSE_STRAND;
    else
        core_.flag &= ~BamRecordImpl::REVERSE_STRAND;
    return *this;
}

BamRecordBuilder& BamRecordBuilder::SetSecondMate(bool ok)
{
    if (ok)
        core_.flag |= BamRecordImpl::MATE_2;
    else
        core_.flag &= ~BamRecordImpl::MATE_2;
    return *this;
}

BamRecordBuilder& BamRecordBuilder::SetSupplementaryAlignment(bool ok)
{
    if (ok)
        core_.flag |= BamRecordImpl::SUPPLEMENTARY;
    else
        core_.flag &= ~BamRecordImpl::SUPPLEMENTARY;
    return *this;
}

}  // namespace BAM
}  // namespace PacBio
