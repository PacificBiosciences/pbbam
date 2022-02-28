#include "PbbamInternalConfig.h"

#include <pbbam/PbiRawData.h>

#include <pbbam/BamFile.h>
#include <pbbam/BamRecord.h>
#include <pbbam/RecordType.h>
#include "PbiIndexIO.h"

#include <boost/numeric/conversion/cast.hpp>

#include <map>
#include <tuple>
#include <type_traits>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace BAM {

// ----------------------------------
// PbiRawBarcodeData implementation
// ----------------------------------

PbiRawBarcodeData::PbiRawBarcodeData(uint32_t numReads)
{
    bcForward_.reserve(numReads);
    bcReverse_.reserve(numReads);
    bcQual_.reserve(numReads);
}

void PbiRawBarcodeData::AddRecord(const BamRecord& b)
{
    // check for any barcode data (both required)
    if (b.HasBarcodes() && b.HasBarcodeQuality()) {

        // fetch data from record
        int16_t bcForward;
        int16_t bcReverse;
        std::tie(bcForward, bcReverse) = b.Barcodes();

        const auto bcQuality = boost::numeric_cast<int8_t>(b.BarcodeQuality());

        // only store actual data if all values >= 0
        if (bcForward >= 0 && bcReverse >= 0 && bcQuality >= 0) {
            bcForward_.push_back(bcForward);
            bcReverse_.push_back(bcReverse);
            bcQual_.push_back(bcQuality);
            return;
        }
    }

    // if we get here, at least one value is either missing or is -1
    bcForward_.push_back(-1);
    bcReverse_.push_back(-1);
    bcQual_.push_back(-1);
}

// ----------------------------------
// PbiRawMappedData implementation
// ----------------------------------

PbiRawMappedData::PbiRawMappedData(uint32_t numReads)
{
    tId_.reserve(numReads);
    tStart_.reserve(numReads);
    tEnd_.reserve(numReads);
    aStart_.reserve(numReads);
    aEnd_.reserve(numReads);
    revStrand_.reserve(numReads);
    nM_.reserve(numReads);
    nMM_.reserve(numReads);
    mapQV_.reserve(numReads);
    nInsOps_.reserve(numReads);
    nDelOps_.reserve(numReads);
}

void PbiRawMappedData::AddRecord(const BamRecord& b)
{
    tId_.push_back(b.ReferenceId());
    tStart_.push_back(b.ReferenceStart());
    tEnd_.push_back(b.ReferenceEnd());
    aStart_.push_back(b.AlignedStart());
    aEnd_.push_back(b.AlignedEnd());
    revStrand_.push_back((b.AlignedStrand() == Strand::REVERSE ? 1 : 0));
    mapQV_.push_back(b.MapQuality());

    const auto matchesAndMismatches = b.NumMatchesAndMismatches();
    nM_.push_back(matchesAndMismatches.first);
    nMM_.push_back(matchesAndMismatches.second);

    const auto insertionAndDeletionOps = b.NumInsertionAndDeletionOperations();
    nInsOps_.push_back(insertionAndDeletionOps.first);
    nDelOps_.push_back(insertionAndDeletionOps.second);
}

uint32_t PbiRawMappedData::NumDeletedBasesAt(size_t recordIndex) const
{
    return NumDeletedAndInsertedBasesAt(recordIndex).first;
}

std::pair<uint32_t, uint32_t> PbiRawMappedData::NumDeletedAndInsertedBasesAt(
    size_t recordIndex) const
{
    const auto aStart = aStart_.at(recordIndex);
    const auto aEnd = aEnd_.at(recordIndex);
    const auto tStart = tStart_.at(recordIndex);
    const auto tEnd = tEnd_.at(recordIndex);
    const auto nM = nM_.at(recordIndex);
    const auto nMM = nMM_.at(recordIndex);

    const auto numIns = (aEnd - aStart - nM - nMM);
    const auto numDel = (tEnd - tStart - nM - nMM);
    return {numDel, numIns};
}

uint32_t PbiRawMappedData::NumInsertedBasesAt(size_t recordIndex) const
{
    return NumDeletedAndInsertedBasesAt(recordIndex).second;
}

// ------------------------------------
// PbiReferenceEntry implementation
// ------------------------------------

const PbiReferenceEntry::ID PbiReferenceEntry::UNMAPPED_ID = static_cast<PbiReferenceEntry::ID>(-1);
const PbiReferenceEntry::Row PbiReferenceEntry::UNSET_ROW = static_cast<PbiReferenceEntry::Row>(-1);

PbiReferenceEntry::PbiReferenceEntry() : PbiReferenceEntry{UNMAPPED_ID, UNSET_ROW, UNSET_ROW} {}

PbiReferenceEntry::PbiReferenceEntry(ID id) : PbiReferenceEntry{id, UNSET_ROW, UNSET_ROW} {}

PbiReferenceEntry::PbiReferenceEntry(ID id, Row beginRow, Row endRow)
    : tId_{id}, beginRow_{beginRow}, endRow_{endRow}
{}

// ------------------------------------
// PbiRawReferenceData implementation
// ------------------------------------

PbiRawReferenceData::PbiRawReferenceData(uint32_t numRefs) { entries_.reserve(numRefs); }

// ----------------------------------
// PbiRawBasicData implementation
// ----------------------------------

PbiRawBasicData::PbiRawBasicData(uint32_t numReads)
{
    rgId_.reserve(numReads);
    qStart_.reserve(numReads);
    qEnd_.reserve(numReads);
    holeNumber_.reserve(numReads);
    readQual_.reserve(numReads);
    ctxtFlag_.reserve(numReads);
    fileOffset_.reserve(numReads);
    fileNumber_.reserve(numReads);
}

void PbiRawBasicData::AddRecord(const BamRecord& b, int64_t offset)
{
    // read group ID
    auto rgId = b.ReadGroupBaseId();
    if (rgId.empty()) {
        rgId = MakeReadGroupId(b.MovieName(), ToString(b.Type()));
    }
    const auto rawid = std::stoul(rgId, nullptr, 16);
    const auto id = static_cast<int32_t>(rawid);
    rgId_.push_back(id);

    // query start/end
    if (IsCcsOrTranscript(b.Type())) {
        qStart_.push_back(0);
        qEnd_.push_back(b.Impl().SequenceLength());
    } else {
        qStart_.push_back(b.QueryStart());
        qEnd_.push_back(b.QueryEnd());
    }

    // add'l basic data
    holeNumber_.push_back(b.HasHoleNumber() ? b.HoleNumber() : 0);
    readQual_.push_back(b.HasReadAccuracy() ? static_cast<float>(b.ReadAccuracy()) : 0.0f);
    ctxtFlag_.push_back(b.HasLocalContextFlags() ? b.LocalContextFlags()
                                                 : Data::LocalContextFlags::NO_LOCAL_CONTEXT);

    // virtual offset of record start
    fileOffset_.push_back(offset);

    // default file number
    fileNumber_.push_back(0);
}

// ----------------------------------
// PbiRawData implementation
// ----------------------------------

PbiRawData::PbiRawData(std::string pbiFilename) : filename_{std::move(pbiFilename)}
{
    PbiIndexIO::LoadFromFile(*this, filename_);
}

PbiRawData::PbiRawData(const DataSet& dataset)
    : sections_{PbiFile::BASIC | PbiFile::MAPPED | PbiFile::BARCODE}
{
    PbiIndexIO::LoadFromDataSet(*this, dataset);
}

const PbiRawBarcodeData& PbiRawData::BarcodeData() const { return barcodeData_; }

PbiRawBarcodeData& PbiRawData::BarcodeData() { return barcodeData_; }

const PbiRawBasicData& PbiRawData::BasicData() const { return basicData_; }

PbiRawBasicData& PbiRawData::BasicData() { return basicData_; }

std::string PbiRawData::Filename() const { return filename_; }

PbiFile::Sections PbiRawData::FileSections() const { return sections_; }

PbiRawData& PbiRawData::FileSections(PbiFile::Sections sections)
{
    sections_ = sections;
    return *this;
}

bool PbiRawData::HasBarcodeData() const { return HasSection(PbiFile::BARCODE); }

bool PbiRawData::HasMappedData() const { return HasSection(PbiFile::MAPPED); }

bool PbiRawData::HasReferenceData() const { return HasSection(PbiFile::REFERENCE); }

bool PbiRawData::HasSection(const PbiFile::Section section) const
{
    return (sections_ & section) != 0;
}

uint32_t PbiRawData::NumReads() const { return numReads_; }

PbiRawData& PbiRawData::NumReads(uint32_t num)
{
    numReads_ = num;
    return *this;
}

const PbiRawMappedData& PbiRawData::MappedData() const { return mappedData_; }

PbiRawMappedData& PbiRawData::MappedData() { return mappedData_; }

const PbiRawReferenceData& PbiRawData::ReferenceData() const { return referenceData_; }

PbiRawReferenceData& PbiRawData::ReferenceData() { return referenceData_; }

PbiFile::VersionEnum PbiRawData::Version() const { return version_; }

PbiRawData& PbiRawData::Version(PbiFile::VersionEnum version)
{
    version_ = version;
    return *this;
}

bool PbiReferenceEntry::operator==(const PbiReferenceEntry& other) const noexcept
{
    return std::tie(tId_, beginRow_, endRow_) ==
           std::tie(other.tId_, other.beginRow_, other.endRow_);
}

// PBI index caching

PbiIndexCache MakePbiIndexCache(const DataSet& dataset)
{
    return MakePbiIndexCache(dataset.BamFiles());
}

PbiIndexCache MakePbiIndexCache(const std::vector<BamFile>& bamFiles)
{
    PbiIndexCache cache = std::make_shared<std::vector<std::shared_ptr<PbiRawData>>>();
    auto& indices = *cache.get();
    for (const auto& bamFile : bamFiles) {
        const auto& pbiFilename = bamFile.PacBioIndexFilename();
        indices.push_back(std::make_shared<PbiRawData>(pbiFilename));
    }
    return cache;
}

PbiIndexCache MakePbiIndexCache(const BamFile& bamFile)
{
    std::vector<BamFile> bamFiles{bamFile};
    return MakePbiIndexCache(bamFiles);
}

}  // namespace BAM
}  // namespace PacBio
