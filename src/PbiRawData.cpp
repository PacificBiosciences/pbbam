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
// Author: Derek Barnett

#include "pbbam/PbiRawData.h"
#include "pbbam/BamFile.h"
#include "pbbam/BamRecord.h"
#include "PbiIndexIO.h"
#include <cassert>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

// ----------------------------------
// PbiRawBarcodeData implementation
// ----------------------------------

PbiRawBarcodeData::PbiRawBarcodeData(void) { }

PbiRawBarcodeData::PbiRawBarcodeData(uint32_t numReads)
{
    bcLeft_.reserve(numReads);
    bcRight_.reserve(numReads);
    bcQual_.reserve(numReads);
    ctxtFlag_.reserve(numReads);
}

PbiRawBarcodeData::PbiRawBarcodeData(const PbiRawBarcodeData& other)
    : bcLeft_(other.bcLeft_)
    , bcRight_(other.bcRight_)
    , bcQual_(other.bcQual_)
    , ctxtFlag_(other.ctxtFlag_)
{ }

PbiRawBarcodeData::PbiRawBarcodeData(PbiRawBarcodeData&& other)
    : bcLeft_(std::move(other.bcLeft_))
    , bcRight_(std::move(other.bcRight_))
    , bcQual_(std::move(other.bcQual_))
    , ctxtFlag_(std::move(other.ctxtFlag_))
{ }

PbiRawBarcodeData& PbiRawBarcodeData::operator=(const PbiRawBarcodeData& other)
{
    bcLeft_ = other.bcLeft_;
    bcRight_ = other.bcRight_;
    bcQual_ = other.bcQual_;
    ctxtFlag_ =other.ctxtFlag_;
    return *this;
}

PbiRawBarcodeData& PbiRawBarcodeData::operator=(PbiRawBarcodeData&& other)
{
    bcLeft_ = std::move(other.bcLeft_);
    bcRight_ = std::move(other.bcRight_);
    bcQual_ = std::move(other.bcQual_);
    ctxtFlag_ = std::move(other.ctxtFlag_);
    return *this;
}

bool PbiRawBarcodeData::AddRecord(const BamRecord& b)
{
    const BamRecordImpl& impl = b.Impl();
    const bool hasBcTag = impl.HasTag("bc");
    const bool hasBqTag = impl.HasTag("bq");
    const bool hasCxTag = impl.HasTag("cx");
    const bool hasBarcodeInfo = hasBcTag && hasBqTag && hasCxTag;
    if (!hasBarcodeInfo)
        return false;

    const vector<uint16_t> bcValue = impl.TagValue("bc").ToUInt16Array();
    assert(bcValue.size() == 2);
    bcLeft_.push_back(bcValue[0]);
    bcRight_.push_back(bcValue[1]);

    bcQual_.push_back(impl.TagValue("bq").ToUInt8());
    ctxtFlag_.push_back(impl.TagValue("cx").ToUInt8());

    return true;
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

bool PbiRawMappedData::AddRecord(const BamRecord& b)
{
    if (!b.IsMapped())
        return false;

    tId_.push_back(b.ReferenceId());
    tStart_.push_back(b.ReferenceStart());
    tEnd_.push_back(b.ReferenceEnd());
    aStart_.push_back(b.AlignedStart());
    aEnd_.push_back(b.AlignedEnd());
    revStrand_.push_back( (b.AlignedStrand() == Strand::REVERSE ? 1 : 0) );
    mapQV_.push_back(b.MapQuality());

//    uint32_t nM = 0;
//    uint32_t nMM = 0;
//    const Cigar& cigar = b.CigarData();
//    auto cigarIter = cigar.cbegin();
//    auto cigarEnd  = cigar.cend();
//    for (; cigarIter != cigarEnd; ++cigarIter) {
//        const CigarOperation& op = (*cigarIter);
//        if (op.Type() == CigarOperationType::SEQUENCE_MATCH)
//            nM += op.Length();
//        else if (op.Type() == CigarOperationType::SEQUENCE_MISMATCH)
//            nMM += op.Length();
//        else if (op.Type() == CigarOperationType::ALIGNMENT_MATCH)
//            throw std::runtime_error("CIGAR operation 'M' is not allowed in PacBio BAM files. Use 'X/=' instead.");
//    }
//    nM_.push_back(nM);
//    nMM_.push_back(nMM);

    nM_.push_back(0);
    nMM_.push_back(0);

    return true;
}
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

PbiRawSubreadData::PbiRawSubreadData(void) { }

PbiRawSubreadData::PbiRawSubreadData(uint32_t numReads)
{
    rgId_.reserve(numReads);
    qStart_.reserve(numReads);
    qEnd_.reserve(numReads);
    holeNumber_.reserve(numReads);
    readQual_.reserve(numReads);
    fileOffset_.reserve(numReads);
}

PbiRawSubreadData::PbiRawSubreadData(const PbiRawSubreadData& other)
    : rgId_(other.rgId_)
    , qStart_(other.qStart_)
    , qEnd_(other.qEnd_)
    , holeNumber_(other.holeNumber_)
    , readQual_(other.readQual_)
    , fileOffset_(other.fileOffset_)
{ }

PbiRawSubreadData::PbiRawSubreadData(PbiRawSubreadData&& other)
    : rgId_(std::move(other.rgId_))
    , qStart_(std::move(other.qStart_))
    , qEnd_(std::move(other.qEnd_))
    , holeNumber_(std::move(other.holeNumber_))
    , readQual_(std::move(other.readQual_))
    , fileOffset_(std::move(other.fileOffset_))
{ }

PbiRawSubreadData& PbiRawSubreadData::operator=(const PbiRawSubreadData& other)
{
    rgId_ = other.rgId_;
    qStart_ = other.qStart_;
    qEnd_ = other.qEnd_;
    holeNumber_ = other.holeNumber_;
    readQual_ = other.readQual_;
    fileOffset_ = other.fileOffset_;
    return *this;
}

PbiRawSubreadData& PbiRawSubreadData::operator=(PbiRawSubreadData&& other)
{
    rgId_ = std::move(other.rgId_);
    qStart_ = std::move(other.qStart_);
    qEnd_ = std::move(other.qEnd_);
    holeNumber_ = std::move(other.holeNumber_);
    readQual_ = std::move(other.readQual_);
    fileOffset_ = std::move(other.fileOffset_);
    return *this;
}

void PbiRawSubreadData::AddRecord(const BamRecord& b, int64_t offset)
{

    string rgId = b.ReadGroupId();
    if (rgId.empty()) {
        // calculate
    }
    const uint32_t rawid = std::stoul(rgId, nullptr, 16);
    const int32_t id = static_cast<int32_t>(rawid);

    rgId_.push_back(id);

    if (b.Type() == RecordType::CCS) {
        qStart_.push_back(-1);
        qEnd_.push_back(-1);
    } else {
        qStart_.push_back(b.QueryStart());
        qEnd_.push_back(b.QueryEnd());
    }

    if (b.HasHoleNumber())
        holeNumber_.push_back(b.HoleNumber());
    else
        holeNumber_.push_back(0); // TODO: what to do?

    if (b.HasReadAccuracy())
        readQual_.push_back(b.ReadAccuracy());
    else
        readQual_.push_back(0); // TODO: what to do?

    fileOffset_.push_back(offset);
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
    : version_(PbiFile::CurrentVersion)
    , sections_(PbiFile::ALL)
    , numReads_(0)
{
    internal::PbiIndexIO::Load(*this, pbiFilename);
}

PbiRawData::PbiRawData(const PbiRawData& other)
    : version_(other.version_)
    , sections_(other.sections_)
    , numReads_(other.numReads_)
    , barcodeData_(other.barcodeData_)
    , mappedData_(other.mappedData_)
    , referenceData_(other.referenceData_)
    , subreadData_(other.subreadData_)
{ }

PbiRawData::PbiRawData(PbiRawData&& other)
    : version_(std::move(other.version_))
    , sections_(std::move(other.sections_))
    , numReads_(std::move(other.numReads_))
    , barcodeData_(std::move(other.barcodeData_))
    , mappedData_(std::move(other.mappedData_))
    , referenceData_(std::move(other.referenceData_))
    , subreadData_(std::move(other.subreadData_))
{ }

PbiRawData& PbiRawData::operator=(const PbiRawData& other)
{
    version_ = other.version_;
    sections_ = other.sections_;
    numReads_ = other.numReads_;
    barcodeData_ = other.barcodeData_;
    mappedData_ = other.mappedData_;
    referenceData_ = other.referenceData_;
    subreadData_ = other.subreadData_;
    return *this;
}

PbiRawData& PbiRawData::operator=(PbiRawData&& other)
{
    version_ = std::move(other.version_);
    sections_ = std::move(other.sections_);
    numReads_ = std::move(other.numReads_);
    barcodeData_ = std::move(other.barcodeData_);
    mappedData_ = std::move(other.mappedData_);
    referenceData_ = std::move(other.referenceData_);
    subreadData_ = std::move(other.subreadData_);
    return *this;
}

PbiRawData::~PbiRawData(void) { }
