// File Description
/// \file PbiRawData.cpp
/// \brief Implements the classes used for working with raw PBI data.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/PbiRawData.h"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <map>

#include <boost/numeric/conversion/cast.hpp>

#include "PbiIndexIO.h"
#include "pbbam/BamFile.h"
#include "pbbam/BamRecord.h"

namespace PacBio {
namespace BAM {
namespace internal {

static std::string ToString(const RecordType type)
{
    // clang-format off
    static const auto lookup = std::map<RecordType, std::string>
    {
        { RecordType::ZMW,        "ZMW" },
        { RecordType::HQREGION,   "HQREGION" },
        { RecordType::SUBREAD,    "SUBREAD" },
        { RecordType::CCS,        "CCS" },
        { RecordType::SCRAP,      "SCRAP" },
        { RecordType::TRANSCRIPT, "TRANSCRIPT" },
        { RecordType::UNKNOWN,    "UNKNOWN" }
    };
    // clang-format on

    try {
        return lookup.at(type);
    } catch (std::exception&) {
        throw std::runtime_error("error: unknown RecordType encountered");
    }
}

}  // namespace internal

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
        const auto barcodes = b.Barcodes();
        const auto barcodeQuality = b.BarcodeQuality();
        const auto bcForward = barcodes.first;
        const auto bcReverse = barcodes.second;
        const auto bcQuality = boost::numeric_cast<int8_t>(barcodeQuality);

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
    return std::make_pair(numDel, numIns);
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

PbiReferenceEntry::PbiReferenceEntry() : tId_(UNMAPPED_ID), beginRow_(UNSET_ROW), endRow_(UNSET_ROW)
{
}

PbiReferenceEntry::PbiReferenceEntry(ID id) : tId_(id), beginRow_(UNSET_ROW), endRow_(UNSET_ROW) {}

PbiReferenceEntry::PbiReferenceEntry(ID id, Row beginRow, Row endRow)
    : tId_(id), beginRow_(beginRow), endRow_(endRow)
{
}

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
    auto rgId = b.ReadGroupId();
    if (rgId.empty()) rgId = MakeReadGroupId(b.MovieName(), internal::ToString(b.Type()));
    const auto rawid = std::stoul(rgId, nullptr, 16);
    const auto id = static_cast<int32_t>(rawid);
    rgId_.push_back(id);

    // query start/end
    if (IsCcsOrTranscript(b.Type())) {
        qStart_.push_back(-1);
        qEnd_.push_back(-1);
    } else {
        qStart_.push_back(b.QueryStart());
        qEnd_.push_back(b.QueryEnd());
    }

    // add'l basic data
    holeNumber_.push_back(b.HasHoleNumber() ? b.HoleNumber() : 0);
    readQual_.push_back(b.HasReadAccuracy() ? static_cast<float>(b.ReadAccuracy()) : 0.0f);
    ctxtFlag_.push_back(b.HasLocalContextFlags() ? b.LocalContextFlags()
                                                 : LocalContextFlags::NO_LOCAL_CONTEXT);

    // virtual offset of record start
    fileOffset_.push_back(offset);

    // default file number
    fileNumber_.push_back(0);
}

// ----------------------------------
// PbiRawData implementation
// ----------------------------------

PbiRawData::PbiRawData(const std::string& pbiFilename) : filename_(pbiFilename)
{
    internal::PbiIndexIO::Load(*this, pbiFilename);
}

PbiRawData::PbiRawData(const DataSet& dataset)
    : sections_(PbiFile::BASIC | PbiFile::MAPPED | PbiFile::BARCODE)
{
    internal::PbiIndexIO::LoadFromDataSet(*this, dataset);
}

}  // namespace BAM
}  // namesapce PacBio
