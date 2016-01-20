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
/// \file PbiBuilder.cpp
/// \brief Implements the PbiBuilder class.
//
// Author: Derek Barnett

#include "pbbam/PbiBuilder.h"
#include "pbbam/BamRecord.h"
#include "pbbam/PbiRawData.h"
#include "MemoryUtils.h"
#include "PbiIndexIO.h"
#include <htslib/bgzf.h>
#include <cassert>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

// -------------------------------------------
// PbiRawReferenceDataBuilder implementation
// -------------------------------------------

// helper for reference data
class PbiRawReferenceDataBuilder
{
public:
    PbiRawReferenceDataBuilder(const size_t numReferenceSequences);

public:
    bool AddRecord(const BamRecord& record,
                   const PbiReferenceEntry::Row rowNumber);
    PbiRawReferenceData Result(void) const;

private:
    int32_t lastRefId_;
    Position lastPos_;
    map<uint32_t, PbiReferenceEntry> rawReferenceEntries_;
};

PbiRawReferenceDataBuilder::PbiRawReferenceDataBuilder(const size_t numReferenceSequences)
    : lastRefId_(-1)
    , lastPos_(-1)
{
    // initialize with number of references we expect to see
    //
    // we can add more later, but want to ensure known references have an entry
    // even if no records are observed mapping to it
    //
    for (size_t i = 0; i < numReferenceSequences; ++i)
        rawReferenceEntries_[i] = PbiReferenceEntry(i);

    // also create an "unmapped" entry
    rawReferenceEntries_[PbiReferenceEntry::UNMAPPED_ID] = PbiReferenceEntry();
}

bool PbiRawReferenceDataBuilder::AddRecord(const BamRecord& record,
                                           const PbiReferenceEntry::Row rowNumber)
{
    // fetch ref ID & pos for record
    const int32_t tId = record.ReferenceId();
    const int32_t pos = record.ReferenceStart();

    // sanity checks to protect against non-coordinate-sorted BAMs
    if (lastRefId_ != tId || (lastRefId_ >= 0 && tId < 0)) {
        if (tId >= 0) {

            // if we've already seen unmapped reads, but our current tId is valid
            //
            // error: unmapped reads should all be at the end (can stop checking refs)
            //
            PbiReferenceEntry& unmappedEntry =
                    rawReferenceEntries_[PbiReferenceEntry::UNMAPPED_ID];
            if (unmappedEntry.beginRow_ != PbiReferenceEntry::UNSET_ROW)
                return false;

            // if we've already seen data for this new tId
            // (remember we're coming from another tId)
            //
            // error: refs are out of order (can stop checking refs)
            //
            PbiReferenceEntry& currentEntry =
                    rawReferenceEntries_[(uint32_t)tId];
            if (currentEntry.beginRow_ != PbiReferenceEntry::UNSET_ROW)
                return false;
        }
        lastRefId_ = tId;
    }
    else if (tId >= 0 && lastPos_ > pos)
        return false; //error: positions out of order

    // update row numbers
    PbiReferenceEntry& entry = rawReferenceEntries_[(uint32_t)tId];
    if (entry.beginRow_ == PbiReferenceEntry::UNSET_ROW)
        entry.beginRow_ = rowNumber;
    entry.endRow_ = rowNumber+1;

    // update pos (for sorting check next go-round)
    lastPos_ = pos;
    return true;
}

PbiRawReferenceData PbiRawReferenceDataBuilder::Result(void) const {
    // PbiReferenceEntries will be sorted thanks to std::map
    // tId will be at end since we're sorting on the uint cast of -1
    PbiRawReferenceData result;
    result.entries_.reserve(rawReferenceEntries_.size());
    auto refIter = rawReferenceEntries_.cbegin();
    const auto refEnd  = rawReferenceEntries_.cend();
    for ( ; refIter != refEnd; ++refIter )
        result.entries_.push_back(refIter->second);
    return result;
}

// ----------------------------------
// PbiBuilderPrivate implementation
// ----------------------------------

class PbiBuilderPrivate
{
public:
    PbiBuilderPrivate(const string& filename,
                      const size_t numReferenceSequences);
    PbiBuilderPrivate(const string& filename,
                      const size_t numReferenceSequences,
                      const bool isCoordinateSorted);
    ~PbiBuilderPrivate(void);

public:
    void AddRecord(const BamRecord& record, const int64_t vOffset);

public:
    bool HasBarcodeData(void) const;
    bool HasMappedData(void) const;
    bool HasReferenceData(void) const;

public:
    unique_ptr<BGZF, HtslibBgzfDeleter> bgzf_;
    PbiRawData rawData_;
    PbiReferenceEntry::Row currentRow_;
    unique_ptr<PbiRawReferenceDataBuilder> refDataBuilder_;
};

PbiBuilderPrivate::PbiBuilderPrivate(const string& filename,
                                     const size_t numReferenceSequences)
    : bgzf_(bgzf_open(filename.c_str(), "wb"))
    , currentRow_(0)
    , refDataBuilder_(nullptr)
{
    if (bgzf_.get() == 0)
        throw std::runtime_error("could not open PBI file for writing");

    if (numReferenceSequences > 0)
        refDataBuilder_.reset(new PbiRawReferenceDataBuilder(numReferenceSequences));
}

PbiBuilderPrivate::PbiBuilderPrivate(const string& filename,
                                     const size_t numReferenceSequences,
                                     const bool isCoordinateSorted)
    : bgzf_(bgzf_open(filename.c_str(), "wb"))
    , currentRow_(0)
    , refDataBuilder_(nullptr)
{
    if (bgzf_.get() == 0)
        throw std::runtime_error("could not open PBI file for writing");

    if (isCoordinateSorted && numReferenceSequences > 0)
        refDataBuilder_.reset(new PbiRawReferenceDataBuilder(numReferenceSequences));
}

PbiBuilderPrivate::~PbiBuilderPrivate(void)
{
    rawData_.NumReads(currentRow_);

    const auto hasBarcodeData   = HasBarcodeData();
    const auto hasMappedData    = HasMappedData();
    const auto hasReferenceData = HasReferenceData();

    // fetch reference data, if available
    if (hasReferenceData) {
        assert(refDataBuilder_);
        rawData_.ReferenceData() = std::move(refDataBuilder_->Result());
    }

    // determine flags
    PbiFile::Sections sections = PbiFile::BASIC;
    if (hasMappedData)    sections |= PbiFile::MAPPED;
    if (hasBarcodeData)   sections |= PbiFile::BARCODE;
    if (hasReferenceData) sections |= PbiFile::REFERENCE;
    rawData_.FileSections(sections);

    // write index contents to file
    BGZF* fp = bgzf_.get();
    PbiIndexIO::WriteHeader(rawData_, fp);
    const uint32_t numReads = rawData_.NumReads();
    if (numReads > 0) {
        PbiIndexIO::WriteBasicData(rawData_.BasicData(), numReads, fp);
        if (hasMappedData)    PbiIndexIO::WriteMappedData(rawData_.MappedData(), numReads, fp);
        if (hasReferenceData) PbiIndexIO::WriteReferenceData(rawData_.ReferenceData(), fp);
        if (hasBarcodeData)   PbiIndexIO::WriteBarcodeData(rawData_.BarcodeData(), numReads, fp);
    }
}

void PbiBuilderPrivate::AddRecord(const BamRecord& record, const int64_t vOffset)
{
    // ensure updated data
    record.ResetCachedPositions();

    // store data
    rawData_.BarcodeData().AddRecord(record);
    rawData_.BasicData().AddRecord(record, vOffset);
    rawData_.MappedData().AddRecord(record);

    if (refDataBuilder_) {

        // stop storing coordinate-sorted reference data if we encounter out-of-order record
        const bool sorted = refDataBuilder_->AddRecord(record, currentRow_);
        if (!sorted)
            refDataBuilder_.reset();
    }

    // increment row counter
    ++currentRow_;
}

bool PbiBuilderPrivate::HasBarcodeData(void) const
{
    // fetch data components
    const auto& barcodeData = rawData_.BarcodeData();
    const auto& bcForward   = barcodeData.bcForward_;
    const auto& bcReverse   = barcodeData.bcReverse_;
    const auto& bcQuality   = barcodeData.bcQual_;

    // ensure valid sizes
    if (bcForward.size() != bcReverse.size() &&
        bcForward.size() != bcQuality.size())
    {
        auto msg = string{ "error: inconsistency in PBI barcode data:\n" };
        msg +=     string{ "  bcForward has " } + to_string(bcForward.size()) + string{ " elements\n" };
        msg +=     string{ "  bcReverse has " } + to_string(bcReverse.size()) + string{ " elements\n" };
        msg +=     string{ "  bcQuality has " } + to_string(bcQuality.size()) + string{ " elements\n" };
        msg +=     string{ "\n" };
        msg +=     string{ "  these containers should contain equal number of elements.\n" };
        throw std::runtime_error(msg);
    }
    assert(bcForward.size() == rawData_.NumReads());

    // check for data
    for (uint32_t i = 0; i < rawData_.NumReads(); ++i) {
        if (bcForward.at(i) != -1 ||
            bcReverse.at(i)  != -1 ||
            bcQuality.at(i)  != -1 )
        {
            return true;
        }
    }
    // no actual data found
    return false;
}

bool PbiBuilderPrivate::HasMappedData(void) const
{
    const auto& mappedData = rawData_.MappedData();
    const auto& tIds = mappedData.tId_;
    assert(tIds.size() == rawData_.NumReads());
    for (const auto tId : tIds) {
        if (tId >= 0)
            return true;
    }
    return false; // all reads unmapped
}

bool PbiBuilderPrivate::HasReferenceData(void) const
{ return bool(refDataBuilder_); }

} // namespace internal
} // namespace BAM
} // namespace PacBio

// ---------------------------
// PbiBuilder implementation
// ---------------------------

PbiBuilder::PbiBuilder(const string& pbiFilename)
    : d_(new internal::PbiBuilderPrivate(pbiFilename, 0))
{ }

PbiBuilder::PbiBuilder(const string& pbiFilename,
                       const size_t numReferenceSequences)
    : d_(new internal::PbiBuilderPrivate(pbiFilename, numReferenceSequences))
{ }

PbiBuilder::PbiBuilder(const string& pbiFilename,
                       const size_t numReferenceSequences,
                       const bool isCoordinateSorted)
    : d_(new internal::PbiBuilderPrivate(pbiFilename,
                                         numReferenceSequences,
                                         isCoordinateSorted))
{ }

PbiBuilder::~PbiBuilder(void) { }

void PbiBuilder::AddRecord(const BamRecord& record, const int64_t vOffset)
{
    internal::BamRecordMemory::UpdateRecordTags(record);
    d_->AddRecord(record, vOffset);
}

const PbiRawData& PbiBuilder::Index(void) const
{ return d_->rawData_; }
