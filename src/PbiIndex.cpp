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
#include "PbiIndex_p.h"
#include "PbiIndexIO.h"

using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {



} // namespace internal
} // namespace BAM
} // namespace PacBio

// ---------------------------
// LookupBase implementation
// ---------------------------

const size_t LookupBase::NullIndex = static_cast<size_t>(-1);

// ----------------------------------
// PerReadLookupBase implementation
// ----------------------------------

// ----------------------------------
// SubreadLookupData implementation
// ----------------------------------

SubreadLookupData::SubreadLookupData(void)
    : PerReadLookupBase()
{ }

SubreadLookupData::SubreadLookupData(const PbiRawSubreadData& rawData)
    : PerReadLookupBase()
    , rgId_(16) // just reserve some space
    , qStart_(MakeLookupMap(rawData.qStart_))
    , qEnd_(MakeLookupMap(rawData.qEnd_))
    , holeNumber_(MakeLookupMap(rawData.holeNumber_))
    , readQual_(MakeLookupMap(rawData.readQual_))
    , fileOffset_(rawData.fileOffset_)
{
    const size_t numElements = rawData.rgId_.size();
    for (size_t i = 0; i < numElements; ++i)
        rgId_[ rawData.rgId_.at(i) ].push_back(i);
}

SubreadLookupData::SubreadLookupData(PbiRawSubreadData&& rawData)
    : PerReadLookupBase()
    , rgId_(16) // just reserve some space
    , qStart_(MakeLookupMap(std::move(rawData.qStart_)))
    , qEnd_(MakeLookupMap(std::move(rawData.qEnd_)))
    , holeNumber_(MakeLookupMap(std::move(rawData.holeNumber_)))
    , readQual_(MakeLookupMap(std::move(rawData.readQual_)))
    , fileOffset_(std::move(rawData.fileOffset_))
{
    const size_t numElements = rawData.rgId_.size();
    for (size_t i = 0; i < numElements; ++i)
        rgId_[ rawData.rgId_.at(i) ].push_back(i);
}

// ----------------------------------
// MappedLookupData implementation
// ----------------------------------

MappedLookupData::MappedLookupData(void)
    : PerReadLookupBase()
{ }

MappedLookupData::MappedLookupData(const PbiRawMappedData& rawData)
    : PerReadLookupBase()
    , tId_(MakeLookupMap(rawData.tId_))
    , tStart_(MakeLookupMap(rawData.tStart_))
    , tEnd_(MakeLookupMap(rawData.tEnd_))
    , aStart_(MakeLookupMap(rawData.aStart_))
    , aEnd_(MakeLookupMap(rawData.aEnd_))
    , nM_(MakeLookupMap(rawData.nM_))
    , nMM_(MakeLookupMap(rawData.nMM_))
    , mapQV_(MakeLookupMap(rawData.mapQV_))
{
    const size_t numElements = rawData.revStrand_.size();
    reverseStrand_.reserve(numElements/2);
    forwardStrand_.reserve(numElements/2);
    for (size_t i = 0; i < numElements; ++i) {
        if (rawData.revStrand_.at(i) == 0)
            forwardStrand_.push_back(i);
        else
            reverseStrand_.push_back(i);
    }
}

MappedLookupData::MappedLookupData(PbiRawMappedData&& rawData)
    : PerReadLookupBase()
    , tId_(MakeLookupMap(std::move(rawData.tId_)))
    , tStart_(MakeLookupMap(std::move(rawData.tStart_)))
    , tEnd_(MakeLookupMap(std::move(rawData.tEnd_)))
    , aStart_(MakeLookupMap(std::move(rawData.aStart_)))
    , aEnd_(MakeLookupMap(std::move(rawData.aEnd_)))
    , nM_(MakeLookupMap(std::move(rawData.nM_)))
    , nMM_(MakeLookupMap(std::move(rawData.nMM_)))
    , mapQV_(MakeLookupMap(std::move(rawData.mapQV_)))
{
    const size_t numElements = rawData.revStrand_.size();
    reverseStrand_.reserve(numElements/2);
    forwardStrand_.reserve(numElements/2);
    for (size_t i = 0; i < numElements; ++i) {
        if (rawData.revStrand_.at(i) == 0)
            forwardStrand_.push_back(i);
        else
            reverseStrand_.push_back(i);
    }
}

// ----------------------------------
// BarcodeLookupData implementation
// ----------------------------------

BarcodeLookupData::BarcodeLookupData(void)
    : PerReadLookupBase()
{ }

BarcodeLookupData::BarcodeLookupData(const PbiRawBarcodeData& rawData)
    : PerReadLookupBase()
    , bcLeft_(MakeLookupMap(rawData.bcLeft_))
    , bcRight_(MakeLookupMap(rawData.bcRight_))
    , bcQual_(MakeLookupMap(rawData.bcQual_))
    , ctxtFlag_(MakeLookupMap(rawData.ctxtFlag_))
{  }

BarcodeLookupData::BarcodeLookupData(PbiRawBarcodeData&& rawData)
    : PerReadLookupBase()
    , bcLeft_(MakeLookupMap(std::move(rawData.bcLeft_)))
    , bcRight_(MakeLookupMap(std::move(rawData.bcRight_)))
    , bcQual_(MakeLookupMap(std::move(rawData.bcQual_)))
    , ctxtFlag_(MakeLookupMap(std::move(rawData.ctxtFlag_)))
{ }

// ----------------------------------
// ReferenceLookupData implementation
// ----------------------------------

ReferenceLookupData::ReferenceLookupData(void)
    : LookupBase()
{ }

ReferenceLookupData::ReferenceLookupData(const PbiRawReferenceData& rawData)
    : LookupBase()
{
    const size_t numEntries = rawData.entries_.size();
    references_.reserve(numEntries);
    for (size_t i = 0; i < numEntries; ++i) {
        const PbiReferenceEntry& entry = rawData.entries_.at(i);
        references_[entry.tId_] = IndexRange(entry.beginRow_, entry.endRow_);
    }
}

ReferenceLookupData::ReferenceLookupData(PbiRawReferenceData&& rawData)
    : LookupBase()
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
{

}

PbiIndexPrivate::PbiIndexPrivate(const PbiRawData& rawIndex)
{

}

PbiIndexPrivate::PbiIndexPrivate(PbiRawData&& rawIndex)
{

}


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

PbiIndex::PbiIndex(const string& pbiFilename)
    : d_(nullptr)
{
    // load raw data from PbiIndexIO
    // update cache

    // not yet implemented
    throw std::exception();
}

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
