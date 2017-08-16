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
#include <cstddef>
#include <cstdint>
#include <map>
#include <cassert>

namespace PacBio {
namespace BAM {
namespace internal {

static
std::string ToString(const RecordType type)
{
    static const auto lookup = std::map<RecordType, std::string>
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

// ------------------------------------
// PbiRawReferenceData implementation
// ------------------------------------

PbiRawReferenceData::PbiRawReferenceData(uint32_t numRefs)
{  entries_.reserve(numRefs); }

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

PbiRawData::PbiRawData(const std::string& pbiFilename)
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

} // namespace BAM
} // namesapce PacBio
