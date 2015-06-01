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

#include "PbiIndexIO.h"
#include "pbbam/BamFile.h"
#include "pbbam/BamRecord.h"
#include "pbbam/EntireFileQuery.h"
#include "MemoryUtils.h"
#include <boost/algorithm/string.hpp>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

// helper for reference data
class PbiRawReferenceDataBuilder
{
public:
    PbiRawReferenceDataBuilder(const BamFile& bam)
        : lastRefId_(-1)
        , lastPos_(-1)
    {
        const BamHeader& header = bam.Header();
        const size_t numReferences = header.Sequences().size();
        for (size_t i = 0; i < numReferences; ++i)
            rawReferenceEntries_[i] = PbiReferenceEntry(i);
        rawReferenceEntries_[PbiReferenceEntry::UNMAPPED_ID] = PbiReferenceEntry();
    }

public:
    bool AddRecord(const BamRecord& record,
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

    PbiRawReferenceData Result(void) const {
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

private:
    int32_t lastRefId_;
    Position lastPos_;
    map<uint32_t, PbiReferenceEntry> rawReferenceEntries_;
};

} // namespace internal
} // namespace BAM
} // namespace PacBio

// ---------------------------
// PbiIndexIO implementation
// ---------------------------

PbiRawData PbiIndexIO::Build(const BamFile& bam)
{
    unique_ptr<samFile,internal::HtslibFileDeleter> htsFile(sam_open(bam.Filename().c_str(), "rb"));
    if (!htsFile)
        throw std::exception();

    unique_ptr<bam_hdr_t, internal::HtslibHeaderDeleter> htsHeader(sam_hdr_read(htsFile.get()));
    if (!htsHeader)
        throw std::exception();

    samFile*   fp  = htsFile.get();
    bam_hdr_t* hdr = htsHeader.get();

    BamRecord record;
    bam1_t* b = internal::BamRecordMemory::GetRawData(record).get();
    if (b == 0)
        throw std::exception();

    // For these optional sections: assume true, we'll mark false if that
    // data type is not present. This allows us to stop checking in during
    // the main loop, and also correctly mark the file at the end.
    bool hasMappedData    = true;
    bool hasBarcodeData   = true;
    bool hasReferenceData = true;

    PbiRawData rawIndex;
    PbiRawBarcodeData& barcodeData = rawIndex.BarcodeData();
    PbiRawMappedData&  mappedData  = rawIndex.MappedData();
    PbiRawSubreadData& subreadData = rawIndex.SubreadData();
    PbiRawReferenceDataBuilder refDataBuilder(bam);

    PbiReferenceEntry::Row rowNumber = 0;
    int64_t offset = bgzf_tell(fp->fp.bgzf);
    while (sam_read1(fp, hdr, b) >= 0) {

        subreadData.AddRecord(record, offset);

        if (hasMappedData)
            hasMappedData &= mappedData.AddRecord(record);

        if (hasReferenceData)
            hasBarcodeData &= barcodeData.AddRecord(record);

        if (hasReferenceData)
            hasReferenceData &= refDataBuilder.AddRecord(record, rowNumber);

        offset = bgzf_tell(fp->fp.bgzf);
        ++rowNumber;
    }
    rawIndex.NumReads(rowNumber);

    // fetch reference data, if available
    if (hasReferenceData)
        rawIndex.ReferenceData() = std::move(refDataBuilder.Result());

    // determine flags
    PbiFile::Sections sections = PbiFile::SUBREAD;
    if (hasMappedData)    sections |= PbiFile::MAPPED;
    if (hasBarcodeData)   sections |= PbiFile::BARCODE;
    if (hasReferenceData) sections |= PbiFile::REFERENCE;
    rawIndex.FileSections(sections);

    return rawIndex;
}

PbiRawData PbiIndexIO::Load(const std::string& pbiFilename)
{
     PbiRawData rawData;
     Load(rawData, pbiFilename);
     return rawData;
}

void PbiIndexIO::Load(PbiRawData& rawData,
                      const string &filename)
{
    // open file for reading
    if (!boost::algorithm::iends_with(filename, ".pbi"))
        throw std::exception();
    std::unique_ptr<BGZF, HtslibBgzfDeleter> bgzf(bgzf_open(filename.c_str(), "rb"));
    BGZF* fp = bgzf.get();
    if (fp == 0)
        throw std::exception();

    // load data
    LoadHeader(rawData, fp);
    const uint32_t numReads = rawData.NumReads();
    if (numReads > 0) {
        LoadSubreadData(rawData.SubreadData(), numReads, fp);
        LoadMappedData(rawData.MappedData(), numReads, fp);
        LoadReferenceData(rawData.ReferenceData(), fp);
        LoadBarcodeData(rawData.BarcodeData(), numReads, fp);
    }
}

void PbiIndexIO::LoadBarcodeData(PbiRawBarcodeData& barcodeData,
                                 const uint32_t numReads,
                                 BGZF* fp)
{
    assert(numReads > 0);

    LoadBgzfVector(fp, barcodeData.bcLeft_,   numReads);
    LoadBgzfVector(fp, barcodeData.bcRight_,  numReads);
    LoadBgzfVector(fp, barcodeData.bcQual_,   numReads);
    LoadBgzfVector(fp, barcodeData.ctxtFlag_, numReads);

    assert(barcodeData.bcLeft_.size()   == numReads);
    assert(barcodeData.bcRight_.size()  == numReads);
    assert(barcodeData.bcQual_.size()   == numReads);
    assert(barcodeData.ctxtFlag_.size() == numReads);
}

void PbiIndexIO::LoadHeader(PbiRawData& index,
                            BGZF* fp)
{
    size_t bytesRead = 0;

    // 'magic' string
    char magic[4];
    bytesRead = bgzf_read(fp, magic, 4);
    if (bytesRead != 4 || strncmp(magic, "PBI\1", 4))
        throw std::exception();

    // version, pbi_flags, & n_reads
    uint32_t version;
    uint16_t sections;
    uint32_t numReads;
    bgzf_read(fp, &version,  sizeof(version));
    bgzf_read(fp, &sections, sizeof(sections));
    bgzf_read(fp, &numReads, sizeof(numReads));
    if (fp->is_be) {
        version  = ed_swap_4(version);
        sections = ed_swap_2(sections);
        numReads = ed_swap_4(numReads);
    }

    index.Version(PbiFile::VersionEnum(version));
    index.FileSections(PbiFile::Sections(sections));
    index.NumReads(numReads);

    // skip reserved section
    size_t reservedLength = 18;
    // adjust depending on version
    char reserved[18];
    bytesRead = bgzf_read(fp, &reserved, reservedLength);
}

void PbiIndexIO::LoadMappedData(PbiRawMappedData& mappedData,
                                const uint32_t numReads,
                                BGZF* fp)
{
    assert(numReads > 0);

    LoadBgzfVector(fp, mappedData.tId_,       numReads);
    LoadBgzfVector(fp, mappedData.tStart_,    numReads);
    LoadBgzfVector(fp, mappedData.tEnd_,      numReads);
    LoadBgzfVector(fp, mappedData.aStart_,    numReads);
    LoadBgzfVector(fp, mappedData.aEnd_,      numReads);
    LoadBgzfVector(fp, mappedData.revStrand_, numReads);
    LoadBgzfVector(fp, mappedData.nM_,        numReads);
    LoadBgzfVector(fp, mappedData.nMM_,       numReads);
    LoadBgzfVector(fp, mappedData.mapQV_,     numReads);

    assert(mappedData.tId_.size()       == numReads);
    assert(mappedData.tStart_.size()    == numReads);
    assert(mappedData.tEnd_.size()      == numReads);
    assert(mappedData.aStart_.size()    == numReads);
    assert(mappedData.aEnd_.size()      == numReads);
    assert(mappedData.revStrand_.size() == numReads);
    assert(mappedData.nM_.size()        == numReads);
    assert(mappedData.nMM_.size()       == numReads);
    assert(mappedData.mapQV_.size()     == numReads);
}

void PbiIndexIO::LoadReferenceData(PbiRawReferenceData& referenceData,
                                   BGZF* fp)
{
    assert(sizeof(PbiReferenceEntry::ID)  == 4);
    assert(sizeof(PbiReferenceEntry::Row) == 4);

    // num refs
    uint32_t numRefs;
    bgzf_read(fp, &numRefs, 4);
    if (fp->is_be)
        numRefs = ed_swap_4(numRefs);

    // reference entries
    referenceData.entries_.clear();
    referenceData.entries_.resize(numRefs);
    for (size_t i = 0; i < numRefs; ++i) {
        PbiReferenceEntry& entry = referenceData.entries_[i];
        bgzf_read(fp, &entry.tId_,      4);
        bgzf_read(fp, &entry.beginRow_, 4);
        bgzf_read(fp, &entry.endRow_,   4);
        if (fp->is_be) {
            entry.tId_      = ed_swap_4(entry.tId_);
            entry.beginRow_ = ed_swap_4(entry.beginRow_);
            entry.endRow_   = ed_swap_4(entry.endRow_);
        }
    }
}

void PbiIndexIO::LoadSubreadData(PbiRawSubreadData& subreadData,
                                 const uint32_t numReads,
                                 BGZF* fp)
{
    assert(numReads > 0);

    LoadBgzfVector(fp, subreadData.rgId_,       numReads);
    LoadBgzfVector(fp, subreadData.qStart_,     numReads);
    LoadBgzfVector(fp, subreadData.qEnd_,       numReads);
    LoadBgzfVector(fp, subreadData.holeNumber_, numReads);
    LoadBgzfVector(fp, subreadData.readQual_,   numReads);
    LoadBgzfVector(fp, subreadData.fileOffset_, numReads);

    assert(subreadData.rgId_.size()       == numReads);
    assert(subreadData.qStart_.size()     == numReads);
    assert(subreadData.qEnd_.size()       == numReads);
    assert(subreadData.holeNumber_.size() == numReads);
    assert(subreadData.readQual_.size()   == numReads);
    assert(subreadData.fileOffset_.size() == numReads);
}

void PbiIndexIO::Save(const PbiRawData& index,
                      const std::string& filename)
{
    std::unique_ptr<BGZF, HtslibBgzfDeleter> bgzf(bgzf_open(filename.c_str(), "wb"));
    BGZF* fp = bgzf.get();
    if (fp == 0)
        throw std::exception();

    WriteHeader(index, fp);
    const uint32_t numReads = index.NumReads();
    if (numReads > 0) {
        WriteSubreadData(index.SubreadData(), numReads, fp);

        if (index.HasMappedData())
            WriteMappedData(index.MappedData(), numReads, fp);
        if (index.HasReferenceData())
            WriteReferenceData(index.ReferenceData(), fp);
        if (index.HasBarcodeData())
            WriteBarcodeData(index.BarcodeData(), numReads, fp);
    }
}

void PbiIndexIO::WriteBarcodeData(const PbiRawBarcodeData& barcodeData,
                                  const uint32_t numReads,
                                  BGZF* fp)
{
    assert(numReads > 0);
    assert(barcodeData.bcLeft_.size()   == numReads);
    assert(barcodeData.bcRight_.size()  == numReads);
    assert(barcodeData.bcQual_.size()   == numReads);
    assert(barcodeData.ctxtFlag_.size() == numReads);

    WriteBgzfVector(fp, barcodeData.bcLeft_);
    WriteBgzfVector(fp, barcodeData.bcRight_);
    WriteBgzfVector(fp, barcodeData.bcQual_);
    WriteBgzfVector(fp, barcodeData.ctxtFlag_);
}

void PbiIndexIO::WriteHeader(const PbiRawData& index,
                             BGZF* fp)
{
    // 'magic' string
    char magic[4];
    strncpy(magic, "PBI\1", 4);
    bgzf_write(fp, magic, 4);

    // version, pbi_flags, & n_reads
    uint32_t version   = static_cast<uint32_t>(index.Version());
    uint16_t pbi_flags = static_cast<uint16_t>(index.FileSections());
    uint32_t numReads  = index.NumReads();
    if (fp->is_be) {
        version   = ed_swap_4(version);
        pbi_flags = ed_swap_2(pbi_flags);
        numReads  = ed_swap_4(numReads);
    }
    bgzf_write(fp, &version,   4);
    bgzf_write(fp, &pbi_flags, 2);
    bgzf_write(fp, &numReads,  4);

    // reserved space
    char reserved[18];
    memset(reserved, 0, 18);
    bgzf_write(fp, reserved, 18);
}

void PbiIndexIO::WriteMappedData(const PbiRawMappedData& mappedData,
                                 const uint32_t numReads,
                                 BGZF* fp)
{
    assert(mappedData.tId_.size()       == numReads);
    assert(mappedData.tStart_.size()    == numReads);
    assert(mappedData.tEnd_.size()      == numReads);
    assert(mappedData.aStart_.size()    == numReads);
    assert(mappedData.aEnd_.size()      == numReads);
    assert(mappedData.revStrand_.size() == numReads);
    assert(mappedData.nM_.size()        == numReads);
    assert(mappedData.nMM_.size()       == numReads);
    assert(mappedData.mapQV_.size()     == numReads);

    WriteBgzfVector(fp, mappedData.tId_);
    WriteBgzfVector(fp, mappedData.tStart_);
    WriteBgzfVector(fp, mappedData.tEnd_);
    WriteBgzfVector(fp, mappedData.aStart_);
    WriteBgzfVector(fp, mappedData.aEnd_);
    WriteBgzfVector(fp, mappedData.revStrand_);
    WriteBgzfVector(fp, mappedData.nM_);
    WriteBgzfVector(fp, mappedData.nMM_);
    WriteBgzfVector(fp, mappedData.mapQV_);
}

void PbiIndexIO::WriteReferenceData(const PbiRawReferenceData& referenceData,
                                    BGZF* fp)
{
    // num_refs
    uint32_t numRefs = referenceData.entries_.size();
    if (fp->is_be)
        numRefs = ed_swap_4(numRefs);
    bgzf_write(fp, &numRefs, 4);

    // reference entries
    numRefs = referenceData.entries_.size(); // need to reset after maybe endian-swapping
    for (size_t i = 0; i < numRefs; ++i) {
        const PbiReferenceEntry& entry = referenceData.entries_[i];
        uint32_t tId      = entry.tId_;
        uint32_t beginRow = entry.beginRow_;
        uint32_t endRow   = entry.endRow_;
        if (fp->is_be) {
            tId      = ed_swap_4(tId);
            beginRow = ed_swap_4(beginRow);
            endRow   = ed_swap_4(endRow);
        }
        bgzf_write(fp, &tId,      4);
        bgzf_write(fp, &beginRow, 4);
        bgzf_write(fp, &endRow,   4);
    }
}

void PbiIndexIO::WriteSubreadData(const PbiRawSubreadData& subreadData,
                                  const uint32_t numReads,
                                  BGZF* fp)
{
    assert(subreadData.rgId_.size()       == numReads);
    assert(subreadData.qStart_.size()     == numReads);
    assert(subreadData.qEnd_.size()       == numReads);
    assert(subreadData.holeNumber_.size() == numReads);
    assert(subreadData.readQual_.size()   == numReads);
    assert(subreadData.fileOffset_.size() == numReads);

    WriteBgzfVector(fp, subreadData.rgId_);
    WriteBgzfVector(fp, subreadData.qStart_);
    WriteBgzfVector(fp, subreadData.qEnd_);
    WriteBgzfVector(fp, subreadData.holeNumber_);
    WriteBgzfVector(fp, subreadData.readQual_);
    WriteBgzfVector(fp, subreadData.fileOffset_);
}
