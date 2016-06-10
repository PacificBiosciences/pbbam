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
/// \file PbiIndex.cpp
/// \brief Implements the PbiIndex class.
//
// Author: Derek Barnett

#include "pbbam/PbiIndex.h"
#include "PbiIndexIO.h"
using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace std;

// ----------------------------------
// SubreadLookupData implementation
// ----------------------------------

BasicLookupData::BasicLookupData(void) { }

BasicLookupData::BasicLookupData(const PbiRawBasicData& rawData)
    : rgId_(rawData.rgId_)
    , qStart_(rawData.qStart_)
    , qEnd_(rawData.qEnd_)
    , holeNumber_(rawData.holeNumber_)
    , readQual_(rawData.readQual_)
    , ctxtFlag_(rawData.ctxtFlag_)
    , fileOffset_(rawData.fileOffset_)
{ }

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

        // nDel, nIns
        const auto indels = rawData.NumDeletedAndInsertedBasesAt(i);
        delRawData[indels.first].push_back(i);
        insRawData[indels.second].push_back(i);

        // strand
        if (rawData.revStrand_.at(i) == 0)
            forwardStrand_.push_back(i);
        else
            reverseStrand_.push_back(i);
    }

    nIns_ = OrderedLookup<uint32_t>(std::move(insRawData));
    nDel_ = OrderedLookup<uint32_t>(std::move(delRawData));
}

// ----------------------------------
// BarcodeLookupData implementation
// ----------------------------------

BarcodeLookupData::BarcodeLookupData(void) { }

BarcodeLookupData::BarcodeLookupData(const PbiRawBarcodeData& rawData)
    : bcForward_(rawData.bcForward_)
    , bcReverse_(rawData.bcReverse_)
    , bcQual_(rawData.bcQual_)

{  }

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

// --------------------------------
// PbiIndexPrivate implementation
// --------------------------------

PbiIndexPrivate::PbiIndexPrivate(void)
    : version_(PbiFile::CurrentVersion)
    , sections_(PbiFile::BASIC)
    , numReads_(0)
{ }

PbiIndexPrivate::PbiIndexPrivate(const PbiRawData& rawIndex)
    : filename_(rawIndex.Filename())
    , version_(rawIndex.Version())
    , sections_(rawIndex.FileSections())
    , numReads_(rawIndex.NumReads())
    , basicData_(rawIndex.BasicData())
    , mappedData_(rawIndex.MappedData())
    , referenceData_(rawIndex.ReferenceData())
    , barcodeData_(rawIndex.BarcodeData())
{ }

PbiIndexPrivate::PbiIndexPrivate(PbiRawData&& rawIndex)
    : filename_(rawIndex.Filename())
    , version_(std::move(rawIndex.Version()))
    , sections_(std::move(rawIndex.FileSections()))
    , numReads_(std::move(rawIndex.NumReads()))
    , basicData_(std::move(rawIndex.BasicData()))
    , mappedData_(std::move(rawIndex.MappedData()))
    , referenceData_(std::move(rawIndex.ReferenceData()))
    , barcodeData_(std::move(rawIndex.BarcodeData()))
{ }

unique_ptr<PbiIndexPrivate> PbiIndexPrivate::DeepCopy(void) const
{
    std::unique_ptr<PbiIndexPrivate> copy(new PbiIndexPrivate);
    copy->filename_ = filename_;
    copy->version_  = version_;
    copy->sections_ = sections_;
    copy->numReads_ = numReads_;
    copy->basicData_     = basicData_;
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
    : d_(std::forward<unique_ptr<PbiIndexPrivate>>(other.d_->DeepCopy()))
{
    // move is ok, since it's a deep-copied, new object
}

PbiIndex::PbiIndex(PbiIndex&& other)
    : d_(std::move(other.d_))
{ }

PbiIndex& PbiIndex::operator=(const PbiIndex& other)
{
    // move is ok, since it's a deep-copied, new object
    d_ = other.d_->DeepCopy();
    return *this;
}

PbiIndex& PbiIndex::operator=(PbiIndex&& other)
{
    d_ = std::move(other.d_);
    return *this;
}

PbiIndex::~PbiIndex(void) { }

string PbiIndex::Filename(void) const
{ return d_->filename_; }
