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

#include "pbbam/PbiIndex.h"
#include "PbiIndexIO.h"
using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace std;

// ----------------------------------
// IndexResultBlock implementation
// ----------------------------------

IndexResultBlock::IndexResultBlock(void)
    : firstIndex_(0)
    , numReads_(0)
    , virtualOffset_(-1)
{ }

IndexResultBlock::IndexResultBlock(size_t idx, size_t numReads)
    : firstIndex_(idx)
    , numReads_(numReads)
    , virtualOffset_(-1)
{ }

bool IndexResultBlock::operator==(const IndexResultBlock& other) const
{
    return firstIndex_ == other.firstIndex_ &&
           numReads_ == other.numReads_ &&
           virtualOffset_ == other.virtualOffset_;
}

bool IndexResultBlock::operator!=(const IndexResultBlock& other) const
{ return !(*this == other); }

// ----------------------------------
// SubreadLookupData implementation
// ----------------------------------

SubreadLookupData::SubreadLookupData(void) { }

SubreadLookupData::SubreadLookupData(const PbiRawSubreadData& rawData)
    : rgId_(rawData.rgId_)
    , qStart_(rawData.qStart_)
    , qEnd_(rawData.qEnd_)
    , holeNumber_(rawData.holeNumber_)
    , readQual_(rawData.readQual_)
    , fileOffset_(rawData.fileOffset_)
{ }

//SubreadLookupData::SubreadLookupData(PbiRawSubreadData&& rawData)
//    : rgId_(std::move(rawData.rgId_))
//    , qStart_(std::move(rawData.qStart_))
//    , qEnd_(std::move(rawData.qEnd_))
//    , holeNumber_(std::move(rawData.holeNumber_))
//    , readQual_(std::move(rawData.readQual_))
//    , fileOffset_(std::move(rawData.fileOffset_))
//{ }

// ----------------------------------
// MappedLookupData implementation
// ----------------------------------

MappedLookupData::MappedLookupData(void) { }

MappedLookupData::MappedLookupData(const PbiRawMappedData& rawData)
    : tId_(rawData.tId_)
    , tStart_(rawData.tStart_)
    , tEnd_(rawData.tEnd_)
    , aStart_(rawData.aStart_)
    , aEnd_(rawData.aEnd_)
    , nM_(rawData.nM_)
    , nMM_(rawData.nMM_)
    , mapQV_(rawData.mapQV_)
{
    const size_t numElements = rawData.revStrand_.size();
    reverseStrand_.reserve(numElements/2);
    forwardStrand_.reserve(numElements/2);

    std::map<uint32_t, IndexList> insRawData;
    std::map<uint32_t, IndexList> delRawData;
    for (size_t i = 0; i < numElements; ++i) {

        // nIns, nDel
        const uint32_t aStart = rawData.aStart_.at(i);
        const uint32_t aEnd   = rawData.aEnd_.at(i);
        const uint32_t tStart = rawData.tStart_.at(i);
        const uint32_t tEnd   = rawData.tEnd_.at(i);
        const uint32_t nM     = rawData.nM_.at(i);
        const uint32_t nMM    = rawData.nMM_.at(i);
        const uint32_t numIns = (aEnd - aStart - nM - nMM);
        const uint32_t numDel = (tEnd - tStart - nM - nMM);
        insRawData[numIns].push_back(i);
        delRawData[numDel].push_back(i);

        // strand
        if (rawData.revStrand_.at(i) == 0)
            forwardStrand_.push_back(i);
        else
            reverseStrand_.push_back(i);
    }

    nIns_ = OrderedLookup<uint32_t>(std::move(insRawData));
    nDel_ = OrderedLookup<uint32_t>(std::move(delRawData));
}

//MappedLookupData::MappedLookupData(PbiRawMappedData&& rawData)
//    : tId_(std::move(rawData.tId_))
//    , tStart_(std::move(rawData.tStart_))
//    , tEnd_(std::move(rawData.tEnd_))
//    , aStart_(std::move(rawData.aStart_))
//    , aEnd_(std::move(rawData.aEnd_))
//    , nM_(std::move(rawData.nM_))
//    , nMM_(std::move(rawData.nMM_))
//    , mapQV_(std::move(rawData.mapQV_))
//{
//    const size_t numElements = rawData.revStrand_.size();
//    reverseStrand_.reserve(numElements/2);
//    forwardStrand_.reserve(numElements/2);
//    for (size_t i = 0; i < numElements; ++i) {
//        if (rawData.revStrand_.at(i) == 0)
//            forwardStrand_.push_back(i);
//        else
//            reverseStrand_.push_back(i);
//    }
//}

// ----------------------------------
// BarcodeLookupData implementation
// ----------------------------------

BarcodeLookupData::BarcodeLookupData(void) { }

BarcodeLookupData::BarcodeLookupData(const PbiRawBarcodeData& rawData)
    : bcLeft_(rawData.bcLeft_)
    , bcRight_(rawData.bcRight_)
    , bcQual_(rawData.bcQual_)
    , ctxtFlag_(rawData.ctxtFlag_)
{  }

//BarcodeLookupData::BarcodeLookupData(PbiRawBarcodeData&& rawData)
//    : bcLeft_(std::move(rawData.bcLeft_))
//    , bcRight_(std::move(rawData.bcRight_))
//    , bcQual_(std::move(rawData.bcQual_))
//    , ctxtFlag_(std::move(rawData.ctxtFlag_))
//{ }

// ----------------------------------
// ReferenceLookupData implementation
// ----------------------------------

ReferenceLookupData::ReferenceLookupData(void) { }

ReferenceLookupData::ReferenceLookupData(const PbiRawReferenceData& rawData)
{
    const size_t numEntries = rawData.entries_.size();
    references_.reserve(numEntries);
    for (size_t i = 0; i < numEntries; ++i) {
        const PbiReferenceEntry& entry = rawData.entries_.at(i);
        references_[entry.tId_] = IndexRange(entry.beginRow_, entry.endRow_);
    }
}

//ReferenceLookupData::ReferenceLookupData(PbiRawReferenceData&& rawData)
//{
//    const size_t numEntries = rawData.entries_.size();
//    references_.reserve(numEntries);
//    for (size_t i = 0; i < numEntries; ++i) {
//        const PbiReferenceEntry& entry = rawData.entries_.at(i);
//        references_[entry.tId_] = IndexRange(entry.beginRow_, entry.endRow_);
//    }
//}

// --------------------------------
// PbiIndexPrivate implementation
// --------------------------------

PbiIndexPrivate::PbiIndexPrivate(void)
    : version_(PbiFile::CurrentVersion)
    , sections_(PbiFile::SUBREAD)
    , numReads_(0)
{ }

PbiIndexPrivate::PbiIndexPrivate(const PbiRawData& rawIndex)
    : version_(rawIndex.Version())
    , sections_(rawIndex.FileSections())
    , numReads_(rawIndex.NumReads())
    , subreadData_(rawIndex.SubreadData())
    , mappedData_(rawIndex.MappedData())
    , referenceData_(rawIndex.ReferenceData())
    , barcodeData_(rawIndex.BarcodeData())
{ }

PbiIndexPrivate::PbiIndexPrivate(PbiRawData&& rawIndex)
    : version_(std::move(rawIndex.Version()))
    , sections_(std::move(rawIndex.FileSections()))
    , numReads_(std::move(rawIndex.NumReads()))
    , subreadData_(std::move(rawIndex.SubreadData()))
    , mappedData_(std::move(rawIndex.MappedData()))
    , referenceData_(std::move(rawIndex.ReferenceData()))
    , barcodeData_(std::move(rawIndex.BarcodeData()))
{ }

unique_ptr<PbiIndexPrivate> PbiIndexPrivate::DeepCopy(void) const
{
    std::unique_ptr<PbiIndexPrivate> copy(new PbiIndexPrivate);
    copy->version_  = version_;
    copy->sections_ = sections_;
    copy->numReads_ = numReads_;
    copy->subreadData_   = subreadData_;
    copy->mappedData_    = mappedData_;
    copy->referenceData_ = referenceData_;
    copy->barcodeData_   = barcodeData_;
    return copy;
}

// -------------------------
// PbiIndex implementation
// -------------------------

PbiIndex::PbiIndex(void)
    : d_(new PbiIndexPrivate)
{ }

PbiIndex::PbiIndex(const string& pbiFilename)
    : d_(new PbiIndexPrivate(PbiRawData(pbiFilename)))
{ }

PbiIndex::PbiIndex(const PbiIndex& other)
    : d_(std::move(other.d_->DeepCopy()))
{
    // move is ok, since it's a deep-copied, new object
}

PbiIndex::PbiIndex(PbiIndex&& other)
    : d_(std::move(other.d_))
{ }

PbiIndex& PbiIndex::operator=(const PbiIndex& other)
{
    // move is ok, since it's a deep-copied, new object
    d_ = std::move(other.d_->DeepCopy());
    return *this;
}

PbiIndex& PbiIndex::operator=(PbiIndex&& other)
{
    d_ = std::move(other.d_);
    return *this;
}

PbiIndex::~PbiIndex(void) { }

PbiFile::Sections PbiIndex::FileSections(void) const
{ return d_->sections_; }

bool PbiIndex::HasBarcodeData(void) const
{ return d_->HasSection(PbiFile::BARCODE); }

bool PbiIndex::HasMappedData(void) const
{ return d_->HasSection(PbiFile::MAPPED); }

bool PbiIndex::HasReferenceData(void) const
{ return d_->HasSection(PbiFile::REFERENCE); }

bool PbiIndex::HasSection(const PbiFile::Section section) const
{ return d_->HasSection(section); }

uint32_t PbiIndex::NumReads(void) const
{ return d_->numReads_; }

PbiFile::VersionEnum PbiIndex::Version(void) const
{ return d_->version_; }

const vector<int64_t>& PbiIndex::VirtualFileOffsets(void) const
{ return d_->subreadData_.fileOffset_; }
