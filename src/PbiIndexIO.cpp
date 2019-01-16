// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "PbiIndexIO.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <boost/algorithm/string.hpp>

#include "MemoryUtils.h"
#include "pbbam/BamFile.h"
#include "pbbam/BamRecord.h"
#include "pbbam/EntireFileQuery.h"
#include "pbbam/MoveAppend.h"
#include "pbbam/PbiBuilder.h"

namespace PacBio {
namespace BAM {
namespace {

// clang-format off
const std::set<PbiFile::Field> AllFields {

    // basic
    PbiFile::Field::RG_ID,
    PbiFile::Field::Q_START,
    PbiFile::Field::Q_END,
    PbiFile::Field::ZMW,
    PbiFile::Field::READ_QUALITY,
    PbiFile::Field::CONTEXT_FLAG,
    PbiFile::Field::VIRTUAL_OFFSET,

    // mapped
    PbiFile::Field::T_ID,
    PbiFile::Field::T_START,
    PbiFile::Field::T_END,
    PbiFile::Field::A_START,
    PbiFile::Field::A_END,
    PbiFile::Field::N_M,
    PbiFile::Field::N_MM,
    PbiFile::Field::N_INS,
    PbiFile::Field::N_DEL,
    PbiFile::Field::MAP_QUALITY,
    PbiFile::Field::STRAND,

    // barcode
    PbiFile::Field::BC_FORWARD,
    PbiFile::Field::BC_REVERSE,
    PbiFile::Field::BC_QUALITY
};
// clang-format on

bool HasField(const std::set<PbiFile::Field>& fields, const PbiFile::Field field)
{
    return fields.find(field) != fields.cend();
}

bool HasSection(const PbiFile::Sections sections, const PbiFile::Section section)
{
    return (sections & section) != 0;
}

}  // anonmyous

PbiIndexIO::PbiIndexIO(const std::string& pbiFilename)
    : pbiFilename_{pbiFilename}, fields_{AllFields}
{
    Open(pbiFilename);
    LoadHeader();
}

PbiIndexIO::PbiIndexIO(const std::string& pbiFilename, const std::set<PbiFile::Field>& fields)
    : pbiFilename_{pbiFilename}, fields_{fields}
{
    Open(pbiFilename);
    LoadHeader();
}

const PbiHeader& PbiIndexIO::Header() const { return header_; }

PbiRawData PbiIndexIO::Load()
{
    // required meta data
    PbiRawData rawData;
    rawData.FileSections(header_.sections);
    rawData.NumReads(header_.numReads);
    rawData.Version(header_.version);
    rawData.filename_ = pbiFilename_;
    if (header_.numReads == 0) return rawData;

    // ensure rewind
    const auto ret = bgzf_seek(fp_.get(), header_.firstRecordOffset, SEEK_SET);
    if (ret != 0) throw std::runtime_error{"could not seek in file: " + pbiFilename_};

    // load from PBI sections
    LoadBasicData(rawData);
    if (HasSection(header_.sections, PbiFile::MAPPED)) LoadMappedData(rawData);
    if (HasSection(header_.sections, PbiFile::REFERENCE)) LoadReferenceData(rawData);
    if (HasSection(header_.sections, PbiFile::BARCODE)) LoadBarcodeData(rawData);
    return rawData;
}

void PbiIndexIO::LoadBarcodeData(PbiRawData& data)
{
    using Field = PbiFile::Field;

    auto& barcodeData = data.BarcodeData();
    MaybeSaveField(barcodeData.bcForward_, Field::BC_FORWARD);
    MaybeSaveField(barcodeData.bcReverse_, Field::BC_REVERSE);
    MaybeSaveField(barcodeData.bcQual_, Field::BC_QUALITY);
}

void PbiIndexIO::LoadBasicData(PbiRawData& data)
{
    using Field = PbiFile::Field;

    auto& basicData = data.BasicData();
    MaybeSaveField(basicData.rgId_, Field::RG_ID);
    MaybeSaveField(basicData.qStart_, Field::Q_START);
    MaybeSaveField(basicData.qEnd_, Field::Q_END);
    MaybeSaveField(basicData.holeNumber_, Field::ZMW);
    MaybeSaveField(basicData.readQual_, Field::READ_QUALITY);
    MaybeSaveField(basicData.ctxtFlag_, Field::CONTEXT_FLAG);
    SaveField(basicData.fileOffset_);  // always store offsets
}

PbiRawData PbiIndexIO::LoadFromDataSet(const DataSet& dataset)
{
    PbiRawData aggregateData;
    AggregateDataSet(aggregateData, dataset);
    return aggregateData;
}

void PbiIndexIO::LoadHeader()
{
    BGZF* bgzf = fp_.get();

    // 'magic' string
    char magic[4];
    auto bytesRead = bgzf_read(bgzf, magic, 4);
    if (bytesRead != 4 || strncmp(magic, "PBI\1", 4))
        throw std::runtime_error{"expected PBI file, found unknown format instead"};

    // read header metadata
    uint32_t version;
    uint16_t sections;
    uint32_t numReads;
    bytesRead = bgzf_read(bgzf, &version, sizeof(version));
    bytesRead = bgzf_read(bgzf, &sections, sizeof(sections));
    bytesRead = bgzf_read(bgzf, &numReads, sizeof(numReads));
    if (bgzf->is_be) {
        version = ed_swap_4(version);
        sections = ed_swap_2(sections);
        numReads = ed_swap_4(numReads);
    }

    // skip reserved space
    constexpr const size_t ReservedLength = 18;
    char reserved[ReservedLength];
    bytesRead = bgzf_read(bgzf, &reserved, ReservedLength);

    // store header info
    header_.version = PbiFile::VersionEnum(version);
    header_.sections = sections;
    header_.numReads = numReads;
    header_.firstRecordOffset = bgzf_tell(fp_);

    // setup dummy buffer
    temp_.resize(numReads);
}

void PbiIndexIO::LoadMappedData(PbiRawData& data)
{
    using Field = PbiFile::Field;

    auto& mappedData = data.MappedData();
    MaybeSaveField(mappedData.tId_, Field::T_ID);
    MaybeSaveField(mappedData.tStart_, Field::T_START);
    MaybeSaveField(mappedData.tEnd_, Field::T_END);
    MaybeSaveField(mappedData.aStart_, Field::A_START);
    MaybeSaveField(mappedData.aEnd_, Field::A_END);
    MaybeSaveField(mappedData.revStrand_, Field::STRAND);
    MaybeSaveField(mappedData.nM_, Field::N_M);
    MaybeSaveField(mappedData.nMM_, Field::N_MM);
    MaybeSaveField(mappedData.mapQV_, Field::MAP_QUALITY);
}

void PbiIndexIO::LoadReferenceData(PbiRawData& data)
{
    assert(sizeof(PbiReferenceEntry::ID) == 4);
    assert(sizeof(PbiReferenceEntry::Row) == 4);

    BGZF* bgzf = fp_.get();

    // num refs
    uint32_t numRefs;
    auto ret = bgzf_read(bgzf, &numRefs, 4);
    if (bgzf->is_be) numRefs = ed_swap_4(numRefs);

    // reference entries
    auto& referenceData = data.ReferenceData();
    referenceData.entries_.clear();
    referenceData.entries_.resize(numRefs);
    for (auto& entry : referenceData.entries_) {
        ret = bgzf_read(bgzf, &entry.tId_, 4);
        ret = bgzf_read(bgzf, &entry.beginRow_, 4);
        ret = bgzf_read(bgzf, &entry.endRow_, 4);
        if (fp_->is_be) {
            entry.tId_ = ed_swap_4(entry.tId_);
            entry.beginRow_ = ed_swap_4(entry.beginRow_);
            entry.endRow_ = ed_swap_4(entry.endRow_);
        }
    }
    UNUSED(ret);
}

template <typename T>
void PbiIndexIO::MaybeSaveField(std::vector<T>& dst, const PbiFile::Field field)
{
    if (HasField(fields_, field)) {
        SaveField(dst);
    } else {
        constexpr const size_t ElementSize = sizeof(T);
        SkipField<ElementSize>();
    }
}

template <typename T>
void PbiIndexIO::SaveField(std::vector<T>& dst)
{
    LoadBgzfVector(fp_.get(), dst, header_.numReads, true);
}

template <size_t ElementSize>
void PbiIndexIO::SkipField()
{
    assert(fp_);
    temp_.resize(header_.numReads);
    auto ret = bgzf_read(fp_.get(), &temp_[0], (header_.numReads * ElementSize));
    UNUSED(ret);
}

void PbiIndexIO::Open(const std::string& filename)
{
    // open file handle
    if (!boost::algorithm::iends_with(filename, ".pbi"))
        throw std::runtime_error{"unsupported file extension on " + filename};
    fp_ = std::unique_ptr<BGZF, HtslibBgzfDeleter>(bgzf_open(filename.c_str(), "rb"));
    if (!fp_) throw std::runtime_error{"could not open PBI file: " + filename + "for reading"};
}

void PbiIndexIO::AggregateDataSet(PbiRawData& aggregateData, const DataSet& dataset)
{
    aggregateData.NumReads(0);
    aggregateData.FileSections(PbiFile::BASIC | PbiFile::MAPPED | PbiFile::BARCODE);
    aggregateData.Version(PbiFile::CurrentVersion);

    const auto bamFiles = dataset.BamFiles();
    uint16_t fileNumber = 0;
    for (const auto& bamFile : bamFiles) {
        PbiRawData currentPbi{bamFile.PacBioIndexFilename()};
        const auto currentPbiCount = currentPbi.NumReads();

        // read count
        aggregateData.NumReads(aggregateData.NumReads() + currentPbiCount);

        // BasicData
        auto& aggregateBasicData = aggregateData.BasicData();
        auto& currentBasicData = currentPbi.BasicData();
        MoveAppend(std::move(currentBasicData.rgId_), aggregateBasicData.rgId_);
        MoveAppend(std::move(currentBasicData.qStart_), aggregateBasicData.qStart_);
        MoveAppend(std::move(currentBasicData.qEnd_), aggregateBasicData.qEnd_);
        MoveAppend(std::move(currentBasicData.holeNumber_), aggregateBasicData.holeNumber_);
        MoveAppend(std::move(currentBasicData.readQual_), aggregateBasicData.readQual_);
        MoveAppend(std::move(currentBasicData.ctxtFlag_), aggregateBasicData.ctxtFlag_);
        MoveAppend(std::move(currentBasicData.fileOffset_), aggregateBasicData.fileOffset_);
        MoveAppend(std::vector<uint16_t>(currentPbiCount, fileNumber),
                   aggregateBasicData.fileNumber_);

        // BarcodeData
        auto& aggregateBarcodeData = aggregateData.BarcodeData();
        if (currentPbi.HasBarcodeData()) {
            auto& currentBarcodeData = currentPbi.BarcodeData();
            MoveAppend(std::move(currentBarcodeData.bcForward_), aggregateBarcodeData.bcForward_);
            MoveAppend(std::move(currentBarcodeData.bcReverse_), aggregateBarcodeData.bcReverse_);
            MoveAppend(std::move(currentBarcodeData.bcQual_), aggregateBarcodeData.bcQual_);
        } else {
            MoveAppend(std::vector<int16_t>(currentPbiCount, -1), aggregateBarcodeData.bcForward_);
            MoveAppend(std::vector<int16_t>(currentPbiCount, -1), aggregateBarcodeData.bcReverse_);
            MoveAppend(std::vector<int8_t>(currentPbiCount, -1), aggregateBarcodeData.bcQual_);
        }

        // MappedData
        auto& aggregateMappedData = aggregateData.MappedData();
        if (currentPbi.HasMappedData()) {
            auto& currentMappedData = currentPbi.MappedData();
            MoveAppend(std::move(currentMappedData.tId_), aggregateMappedData.tId_);
            MoveAppend(std::move(currentMappedData.tStart_), aggregateMappedData.tStart_);
            MoveAppend(std::move(currentMappedData.tEnd_), aggregateMappedData.tEnd_);
            MoveAppend(std::move(currentMappedData.aStart_), aggregateMappedData.aStart_);
            MoveAppend(std::move(currentMappedData.aEnd_), aggregateMappedData.aEnd_);
            MoveAppend(std::move(currentMappedData.revStrand_), aggregateMappedData.revStrand_);
            MoveAppend(std::move(currentMappedData.nM_), aggregateMappedData.nM_);
            MoveAppend(std::move(currentMappedData.nMM_), aggregateMappedData.nMM_);
            MoveAppend(std::move(currentMappedData.mapQV_), aggregateMappedData.mapQV_);
        } else {
            MoveAppend(std::vector<int32_t>(currentPbiCount, -1), aggregateMappedData.tId_);
            MoveAppend(std::vector<uint32_t>(currentPbiCount, UnmappedPosition),
                       aggregateMappedData.tStart_);
            MoveAppend(std::vector<uint32_t>(currentPbiCount, UnmappedPosition),
                       aggregateMappedData.tEnd_);
            MoveAppend(std::vector<uint32_t>(currentPbiCount, UnmappedPosition),
                       aggregateMappedData.aStart_);
            MoveAppend(std::vector<uint32_t>(currentPbiCount, UnmappedPosition),
                       aggregateMappedData.aEnd_);
            MoveAppend(std::vector<uint8_t>(currentPbiCount, 0), aggregateMappedData.revStrand_);
            MoveAppend(std::vector<uint32_t>(currentPbiCount, 0), aggregateMappedData.nM_);
            MoveAppend(std::vector<uint32_t>(currentPbiCount, 0), aggregateMappedData.nMM_);
            MoveAppend(std::vector<uint8_t>(currentPbiCount, 255), aggregateMappedData.mapQV_);
        }

        ++fileNumber;
    }
}

}  // namespace BAM
}  // namespace PacBio
