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
/// \file PbiRawData.cpp
/// \brief Implements the classes used for working with raw PBI data.
//
// Author: Derek Barnett

#include "pbbam/PbiRawData.h"
#include "pbbam/BamFile.h"
#include "pbbam/BamRecord.h"
#include "PbiIndexIO.h"
#include <boost/numeric/conversion/cast.hpp>
#include <map>
#include <cassert>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

static
string ToString(const RecordType type)
{
    static const auto lookup = map<RecordType, string>
    {
        { RecordType::ZMW,        "ZMW" },
        { RecordType::HQREGION,   "HQREGION" },
        { RecordType::SUBREAD,    "SUBREAD" },
        { RecordType::CCS,        "CCS" },
        { RecordType::SCRAP,      "SCRAP" },
        { RecordType::UNKNOWN,    "UNKNOWN" }
    };

    try {
        return lookup.at(type);
    } catch (std::exception&) {
        throw std::runtime_error("error: unknown RecordType encountered");
    }
}

} // namespace internal
} // namespace BAM
} // namesapce PacBio

// ----------------------------------
// PbiRawBarcodeData implementation
// ----------------------------------

PbiRawBarcodeData::PbiRawBarcodeData(void) { }

PbiRawBarcodeData::PbiRawBarcodeData(uint32_t numReads)
{
    bcForward_.reserve(numReads);
    bcReverse_.reserve(numReads);
    bcQual_.reserve(numReads);
}

PbiRawBarcodeData::PbiRawBarcodeData(const PbiRawBarcodeData& other)
    : bcForward_(other.bcForward_)
    , bcReverse_(other.bcReverse_)
    , bcQual_(other.bcQual_)
{ }

PbiRawBarcodeData::PbiRawBarcodeData(PbiRawBarcodeData&& other)
    : bcForward_(std::move(other.bcForward_))
    , bcReverse_(std::move(other.bcReverse_))
    , bcQual_(std::move(other.bcQual_))
{ }

PbiRawBarcodeData& PbiRawBarcodeData::operator=(const PbiRawBarcodeData& other)
{
    bcForward_ = other.bcForward_;
    bcReverse_ = other.bcReverse_;
    bcQual_ = other.bcQual_;
    return *this;
}

PbiRawBarcodeData& PbiRawBarcodeData::operator=(PbiRawBarcodeData&& other)
{
    bcForward_ = std::move(other.bcForward_);
    bcReverse_ = std::move(other.bcReverse_);
    bcQual_ = std::move(other.bcQual_);
    return *this;
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
        if (bcForward >= 0 && bcReverse >=0 && bcQuality >= 0) {
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

PbiRawMappedData::PbiRawMappedData(void) { }

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

PbiRawMappedData::PbiRawMappedData(const PbiRawMappedData& other)
    : tId_(other.tId_)
    , tStart_(other.tStart_)
    , tEnd_(other.tEnd_)
    , aStart_(other.aStart_)
    , aEnd_(other.aEnd_)
    , revStrand_(other.revStrand_)
    , nM_(other.nM_)
    , nMM_(other.nMM_)
    , mapQV_(other.mapQV_)
{ }

PbiRawMappedData::PbiRawMappedData(PbiRawMappedData&& other)
    : tId_(std::move(other.tId_))
    , tStart_(std::move(other.tStart_))
    , tEnd_(std::move(other.tEnd_))
    , aStart_(std::move(other.aStart_))
    , aEnd_(std::move(other.aEnd_))
    , revStrand_(std::move(other.revStrand_))
    , nM_(std::move(other.nM_))
    , nMM_(std::move(other.nMM_))
    , mapQV_(std::move(other.mapQV_))
{ }

PbiRawMappedData& PbiRawMappedData::operator=(const PbiRawMappedData& other)
{
    tId_ = other.tId_;
    tStart_ = other.tStart_;
    tEnd_ = other.tEnd_;
    aStart_ = other.aStart_;
    aEnd_ = other.aEnd_;
    revStrand_ = other.revStrand_;
    nM_ = other.nM_;
    nMM_ = other.nMM_;
    mapQV_ = other.mapQV_;
    return *this;
}

PbiRawMappedData& PbiRawMappedData::operator=(PbiRawMappedData&& other)
{
    tId_ = std::move(other.tId_);
    tStart_ = std::move(other.tStart_);
    tEnd_ = std::move(other.tEnd_);
    aStart_ = std::move(other.aStart_);
    aEnd_ = std::move(other.aEnd_);
    revStrand_ = std::move(other.revStrand_);
    nM_ = std::move(other.nM_);
    nMM_ = std::move(other.nMM_);
    mapQV_ = std::move(other.mapQV_);
    return *this;
}

void PbiRawMappedData::AddRecord(const BamRecord& b)
{
    tId_.push_back(b.ReferenceId());
    tStart_.push_back(b.ReferenceStart());
    tEnd_.push_back(b.ReferenceEnd());
    aStart_.push_back(b.AlignedStart());
    aEnd_.push_back(b.AlignedEnd());
    revStrand_.push_back( (b.AlignedStrand() == Strand::REVERSE ? 1 : 0) );
    mapQV_.push_back(b.MapQuality());

    const auto matchesAndMismatches = b.NumMatchesAndMismatches();
    nM_.push_back(matchesAndMismatches.first);
    nMM_.push_back(matchesAndMismatches.second);
}

uint32_t PbiRawMappedData::NumDeletedBasesAt(size_t recordIndex) const
{ return NumDeletedAndInsertedBasesAt(recordIndex).first; }

std::pair<uint32_t, uint32_t> PbiRawMappedData::NumDeletedAndInsertedBasesAt(size_t recordIndex) const
{
    const auto aStart = aStart_.at(recordIndex);
    const auto aEnd   = aEnd_.at(recordIndex);
    const auto tStart = tStart_.at(recordIndex);
    const auto tEnd   = tEnd_.at(recordIndex);
    const auto nM     = nM_.at(recordIndex);
    const auto nMM    = nMM_.at(recordIndex);
    const auto numIns = (aEnd - aStart - nM - nMM);
    const auto numDel = (tEnd - tStart - nM - nMM);
    return std::make_pair(numDel, numIns);
}

uint32_t PbiRawMappedData::NumInsertedBasesAt(size_t recordIndex) const
{ return NumDeletedAndInsertedBasesAt(recordIndex).second; }

// ------------------------------------
// PbiReferenceEntry implementation
// ------------------------------------

const PbiReferenceEntry::ID  PbiReferenceEntry::UNMAPPED_ID = static_cast<PbiReferenceEntry::ID>(-1);
const PbiReferenceEntry::Row PbiReferenceEntry::UNSET_ROW   = static_cast<PbiReferenceEntry::Row>(-1);

PbiReferenceEntry::PbiReferenceEntry(void)
    : tId_(UNMAPPED_ID)
    , beginRow_(UNSET_ROW)
    , endRow_(UNSET_ROW)
{ }

PbiReferenceEntry::PbiReferenceEntry(ID id)
    : tId_(id)
    , beginRow_(UNSET_ROW)
    , endRow_(UNSET_ROW)
{ }

PbiReferenceEntry::PbiReferenceEntry(ID id, Row beginRow, Row endRow)
    : tId_(id)
    , beginRow_(beginRow)
    , endRow_(endRow)
{ }

PbiReferenceEntry::PbiReferenceEntry(const PbiReferenceEntry& other)
    : tId_(other.tId_)
    , beginRow_(other.beginRow_)
    , endRow_(other.endRow_)
{ }

PbiReferenceEntry::PbiReferenceEntry(PbiReferenceEntry&& other)
    : tId_(std::move(other.tId_))
    , beginRow_(std::move(other.beginRow_))
    , endRow_(std::move(other.endRow_))
{ }

PbiReferenceEntry& PbiReferenceEntry::operator=(const PbiReferenceEntry& other)
{
    tId_ = other.tId_;
    beginRow_ = other.beginRow_;
    endRow_ = other.endRow_;
    return *this;
}

PbiReferenceEntry& PbiReferenceEntry::operator=(PbiReferenceEntry&& other)
{
    tId_ = std::move(other.tId_);
    beginRow_ = std::move(other.beginRow_);
    endRow_ = std::move(other.endRow_);
    return *this;
}

// ------------------------------------
// PbiRawReferenceData implementation
// ------------------------------------

PbiRawReferenceData::PbiRawReferenceData(void) { }

PbiRawReferenceData::PbiRawReferenceData(uint32_t numRefs)
{  entries_.reserve(numRefs); }

PbiRawReferenceData::PbiRawReferenceData(const PbiRawReferenceData& other)
    : entries_(other.entries_)
{ }

PbiRawReferenceData::PbiRawReferenceData(PbiRawReferenceData&& other)
    : entries_(std::move(other.entries_))
{ }

PbiRawReferenceData& PbiRawReferenceData::operator=(const PbiRawReferenceData& other)
{
    entries_ = other.entries_;
    return *this;
}

PbiRawReferenceData& PbiRawReferenceData::operator=(PbiRawReferenceData&& other)
{
    entries_ = std::move(other.entries_);
    return *this;
}

// ----------------------------------
// PbiRawSubreadData implementation
// ----------------------------------

PbiRawBasicData::PbiRawBasicData(void) { }

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

PbiRawBasicData::PbiRawBasicData(const PbiRawBasicData& other)
    : rgId_(other.rgId_)
    , qStart_(other.qStart_)
    , qEnd_(other.qEnd_)
    , holeNumber_(other.holeNumber_)
    , readQual_(other.readQual_)
    , ctxtFlag_(other.ctxtFlag_)
    , fileOffset_(other.fileOffset_)
    , fileNumber_(other.fileNumber_)
{ }

PbiRawBasicData::PbiRawBasicData(PbiRawBasicData&& other)
    : rgId_(std::move(other.rgId_))
    , qStart_(std::move(other.qStart_))
    , qEnd_(std::move(other.qEnd_))
    , holeNumber_(std::move(other.holeNumber_))
    , readQual_(std::move(other.readQual_))
    , ctxtFlag_(std::move(other.ctxtFlag_))
    , fileOffset_(std::move(other.fileOffset_))
    , fileNumber_(std::move(other.fileNumber_))
{ }

PbiRawBasicData& PbiRawBasicData::operator=(const PbiRawBasicData& other)
{
    rgId_ = other.rgId_;
    qStart_ = other.qStart_;
    qEnd_ = other.qEnd_;
    holeNumber_ = other.holeNumber_;
    readQual_ = other.readQual_;
    ctxtFlag_ = other.ctxtFlag_;
    fileOffset_ = other.fileOffset_;
    fileNumber_ = other.fileNumber_;
    return *this;
}

PbiRawBasicData& PbiRawBasicData::operator=(PbiRawBasicData&& other)
{
    rgId_ = std::move(other.rgId_);
    qStart_ = std::move(other.qStart_);
    qEnd_ = std::move(other.qEnd_);
    holeNumber_ = std::move(other.holeNumber_);
    readQual_ = std::move(other.readQual_);
    ctxtFlag_ = std::move(other.ctxtFlag_);
    fileOffset_ = std::move(other.fileOffset_);
    fileNumber_ = std::move(other.fileNumber_);
    return *this;
}

void PbiRawBasicData::AddRecord(const BamRecord& b, int64_t offset)
{
    // read group ID
    auto rgId = b.ReadGroupId();
    if (rgId.empty())
        rgId = MakeReadGroupId(b.MovieName(), internal::ToString(b.Type()));
    const uint32_t rawid = std::stoul(rgId, nullptr, 16);
    const int32_t id = static_cast<int32_t>(rawid);
    rgId_.push_back(id);

    // query start/end
    if (b.Type() == RecordType::CCS) {
        qStart_.push_back(-1);
        qEnd_.push_back(-1);
    } else {
        qStart_.push_back(b.QueryStart());
        qEnd_.push_back(b.QueryEnd());
    }

    // add'l basic data
    holeNumber_.push_back(b.HasHoleNumber() ? b.HoleNumber() : 0);
    readQual_.push_back(b.HasReadAccuracy() ? static_cast<float>(b.ReadAccuracy()) : 0.0f);
    ctxtFlag_.push_back(b.HasLocalContextFlags() ? b.LocalContextFlags() : LocalContextFlags::NO_LOCAL_CONTEXT);

    // virtual offset of record start
    fileOffset_.push_back(offset);

    // default file number
    fileNumber_.push_back(0);
}

// ----------------------------------
// PbiRawData implementation
// ----------------------------------

PbiRawData::PbiRawData(void)
    : version_(PbiFile::CurrentVersion)
    , sections_(PbiFile::ALL)
    , numReads_(0)
{ }

PbiRawData::PbiRawData(const string& pbiFilename)
    : filename_(pbiFilename)
    , version_(PbiFile::CurrentVersion)
    , sections_(PbiFile::ALL)
    , numReads_(0)
{
    internal::PbiIndexIO::Load(*this, pbiFilename);
}

PbiRawData::PbiRawData(const DataSet& dataset)
    : version_(PbiFile::CurrentVersion)
    , sections_(PbiFile::BASIC | PbiFile::MAPPED | PbiFile::BARCODE)
    , numReads_(0)
{
    internal::PbiIndexIO::LoadFromDataSet(*this, dataset);
}

PbiRawData::PbiRawData(const PbiRawData& other)
    : filename_(other.filename_)
    , version_(other.version_)
    , sections_(other.sections_)
    , numReads_(other.numReads_)
    , barcodeData_(other.barcodeData_)
    , mappedData_(other.mappedData_)
    , referenceData_(other.referenceData_)
    , basicData_(other.basicData_)
{ }

PbiRawData::PbiRawData(PbiRawData&& other)
    : filename_(std::move(other.filename_))
    , version_(std::move(other.version_))
    , sections_(std::move(other.sections_))
    , numReads_(std::move(other.numReads_))
    , barcodeData_(std::move(other.barcodeData_))
    , mappedData_(std::move(other.mappedData_))
    , referenceData_(std::move(other.referenceData_))
    , basicData_(std::move(other.basicData_))
{ }

PbiRawData& PbiRawData::operator=(const PbiRawData& other)
{
    filename_ = other.filename_;
    version_ = other.version_;
    sections_ = other.sections_;
    numReads_ = other.numReads_;
    barcodeData_ = other.barcodeData_;
    mappedData_ = other.mappedData_;
    referenceData_ = other.referenceData_;
    basicData_ = other.basicData_;
    return *this;
}

PbiRawData& PbiRawData::operator=(PbiRawData&& other)
{
    filename_ = std::move(other.filename_);
    version_ = std::move(other.version_);
    sections_ = std::move(other.sections_);
    numReads_ = std::move(other.numReads_);
    barcodeData_ = std::move(other.barcodeData_);
    mappedData_ = std::move(other.mappedData_);
    referenceData_ = std::move(other.referenceData_);
    basicData_ = std::move(other.basicData_);
    return *this;
}

PbiRawData::~PbiRawData(void) { }
