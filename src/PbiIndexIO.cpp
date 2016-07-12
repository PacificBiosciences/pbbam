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
#include "pbbam/PbiBuilder.h"
#include "MemoryUtils.h"
#include <boost/algorithm/string.hpp>
#include <vector>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

// \brief Appends content of src vector to dst vector using move semantics.
///
/// \param[in]     src  Input vector that will be empty after execution
/// \param[in,out] dst  Output vector that will be appended to
///
template <typename T>
inline void MoveAppend(std::vector<T>& src, std::vector<T>& dst) noexcept
{
    if (dst.empty())
    {
        dst = std::move(src);
    }
    else
    {
        dst.reserve(dst.size() + src.size());
        std::move(src.begin(), src.end(), std::back_inserter(dst));
        src.clear();
    }
}

/// \brief Appends content of src vector to dst vector using move semantics.
///
/// \param[in]     src  Input vector via perfect forwarding
/// \param[in,out] dst  Output vector that will be appended to
///
template <typename T>
inline void MoveAppend(std::vector<T>&& src, std::vector<T>& dst) noexcept
{
    if (dst.empty())
    {
        dst = std::move(src);
    }
    else
    {
        dst.reserve(dst.size() + src.size());
        std::move(src.begin(), src.end(), std::back_inserter(dst));
        src.clear();
    }
}

} // namespace internal
} // namespace BAM
} // namespace PacBio

// ---------------------------
// PbiIndexIO implementation
// ---------------------------

PbiRawData PbiIndexIO::Load(const std::string& pbiFilename)
{
     PbiRawData rawData;
     Load(rawData, pbiFilename);
     return rawData;
}

void PbiIndexIO::Load(PbiRawData& rawData,
                      const string& filename)
{
    // open file for reading
    if (!boost::algorithm::iends_with(filename, ".pbi"))
        throw std::runtime_error("unsupported file extension");
    std::unique_ptr<BGZF, HtslibBgzfDeleter> bgzf(bgzf_open(filename.c_str(), "rb"));
    BGZF* fp = bgzf.get();
    if (fp == 0)
        throw std::runtime_error("could not open PBI file for reading");

    // load data
    LoadHeader(rawData, fp);
    const uint32_t numReads = rawData.NumReads();
    if (numReads > 0) {
        LoadBasicData(rawData.BasicData(), numReads, fp);
        if (rawData.HasMappedData())
            LoadMappedData(rawData.MappedData(), numReads, fp);
        if (rawData.HasReferenceData())
            LoadReferenceData(rawData.ReferenceData(), fp);
        if (rawData.HasBarcodeData())
            LoadBarcodeData(rawData.BarcodeData(), numReads, fp);
    }
}

void PbiIndexIO::LoadFromDataSet(PbiRawData& aggregateData,
                                 const DataSet& dataset)
{
    aggregateData.NumReads(0);
    aggregateData.FileSections(PbiFile::BASIC | PbiFile::MAPPED | PbiFile::BARCODE);
    aggregateData.Version(PbiFile::CurrentVersion);

    const auto bamFiles = dataset.BamFiles();
    uint16_t fileNumber = 0;
    for (const auto& bamFile : bamFiles) {
        PbiRawData currentPbi{bamFile.PacBioIndexFilename()};
        const size_t currentPbiCount = currentPbi.NumReads();

        // read count
        aggregateData.NumReads(aggregateData.NumReads()+currentPbiCount);

        // BasicData
        PbiRawBasicData& aggregateBasicData = aggregateData.BasicData();
        PbiRawBasicData& currentBasicData   = currentPbi.BasicData();
        MoveAppend(std::move(currentBasicData.rgId_),       aggregateBasicData.rgId_);
        MoveAppend(std::move(currentBasicData.qStart_),     aggregateBasicData.qStart_);
        MoveAppend(std::move(currentBasicData.qEnd_),       aggregateBasicData.qEnd_);
        MoveAppend(std::move(currentBasicData.holeNumber_), aggregateBasicData.holeNumber_);
        MoveAppend(std::move(currentBasicData.readQual_),   aggregateBasicData.readQual_);
        MoveAppend(std::move(currentBasicData.ctxtFlag_),   aggregateBasicData.ctxtFlag_);
        MoveAppend(std::move(currentBasicData.fileOffset_), aggregateBasicData.fileOffset_);
        MoveAppend(std::vector<uint16_t>(currentPbiCount, fileNumber), aggregateBasicData.fileNumber_);

        // BarcodeData
        PbiRawBarcodeData& aggregateBarcodeData = aggregateData.BarcodeData();
        if (currentPbi.HasBarcodeData()) {
            PbiRawBarcodeData& currentBarcodeData  = currentPbi.BarcodeData();
            MoveAppend(std::move(currentBarcodeData.bcForward_), aggregateBarcodeData.bcForward_);
            MoveAppend(std::move(currentBarcodeData.bcReverse_), aggregateBarcodeData.bcReverse_);
            MoveAppend(std::move(currentBarcodeData.bcQual_),    aggregateBarcodeData.bcQual_);
        } else {
            MoveAppend(std::vector<int16_t>(currentPbiCount, -1), aggregateBarcodeData.bcForward_);
            MoveAppend(std::vector<int16_t>(currentPbiCount, -1), aggregateBarcodeData.bcReverse_);
            MoveAppend(std::vector<int8_t>(currentPbiCount, -1),  aggregateBarcodeData.bcQual_);
        }

        // MappedData
        PbiRawMappedData& aggregateMappedData = aggregateData.MappedData();
        if (currentPbi.HasMappedData()) {
            PbiRawMappedData& currentMappedData  = currentPbi.MappedData();
            MoveAppend(std::move(currentMappedData.tId_),       aggregateMappedData.tId_);
            MoveAppend(std::move(currentMappedData.tStart_),    aggregateMappedData.tStart_);
            MoveAppend(std::move(currentMappedData.tEnd_),      aggregateMappedData.tEnd_);
            MoveAppend(std::move(currentMappedData.aStart_),    aggregateMappedData.aStart_);
            MoveAppend(std::move(currentMappedData.aEnd_),      aggregateMappedData.aEnd_);
            MoveAppend(std::move(currentMappedData.revStrand_), aggregateMappedData.revStrand_);
            MoveAppend(std::move(currentMappedData.nM_),        aggregateMappedData.nM_);
            MoveAppend(std::move(currentMappedData.nMM_),       aggregateMappedData.nMM_);
            MoveAppend(std::move(currentMappedData.mapQV_),     aggregateMappedData.mapQV_);
        } else {
            MoveAppend(std::vector<int32_t>(currentPbiCount, -1), aggregateMappedData.tId_);
            MoveAppend(std::vector<uint32_t>(currentPbiCount, UnmappedPosition), aggregateMappedData.tStart_);
            MoveAppend(std::vector<uint32_t>(currentPbiCount, UnmappedPosition), aggregateMappedData.tEnd_);
            MoveAppend(std::vector<uint32_t>(currentPbiCount, UnmappedPosition), aggregateMappedData.aStart_);
            MoveAppend(std::vector<uint32_t>(currentPbiCount, UnmappedPosition), aggregateMappedData.aEnd_);
            MoveAppend(std::vector<uint8_t>(currentPbiCount, 0),   aggregateMappedData.revStrand_);
            MoveAppend(std::vector<uint32_t>(currentPbiCount, 0),  aggregateMappedData.nM_);
            MoveAppend(std::vector<uint32_t>(currentPbiCount, 0),  aggregateMappedData.nMM_);
            MoveAppend(std::vector<uint8_t>(currentPbiCount, 255), aggregateMappedData.mapQV_);
        }

        ++fileNumber;
    }
}

void PbiIndexIO::LoadBarcodeData(PbiRawBarcodeData& barcodeData,
                                 const uint32_t numReads,
                                 BGZF* fp)
{
    assert(numReads > 0);
    (void)numReads; // quash warnings building in release mode

    LoadBgzfVector(fp, barcodeData.bcForward_, numReads);
    LoadBgzfVector(fp, barcodeData.bcReverse_, numReads);
    LoadBgzfVector(fp, barcodeData.bcQual_,    numReads);

    assert(barcodeData.bcForward_.size() == numReads);
    assert(barcodeData.bcReverse_.size() == numReads);
    assert(barcodeData.bcQual_.size()    == numReads);
}

void PbiIndexIO::LoadHeader(PbiRawData& index,
                            BGZF* fp)
{
    size_t bytesRead = 0;

    // 'magic' string
    char magic[4];
    bytesRead = bgzf_read(fp, magic, 4);
    if (bytesRead != 4 || strncmp(magic, "PBI\1", 4))
        throw std::runtime_error("expected PBI file, found unknown format instead");

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
    (void)numReads; // quash warnings building in release mode

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

void PbiIndexIO::LoadBasicData(PbiRawBasicData& basicData,
                                 const uint32_t numReads,
                                 BGZF* fp)
{
    assert(numReads > 0);
    (void)numReads; // quash warnings building in release mode

    LoadBgzfVector(fp, basicData.rgId_,       numReads);
    LoadBgzfVector(fp, basicData.qStart_,     numReads);
    LoadBgzfVector(fp, basicData.qEnd_,       numReads);
    LoadBgzfVector(fp, basicData.holeNumber_, numReads);
    LoadBgzfVector(fp, basicData.readQual_,   numReads);
    LoadBgzfVector(fp, basicData.ctxtFlag_,   numReads);
    LoadBgzfVector(fp, basicData.fileOffset_, numReads);

    assert(basicData.rgId_.size()       == numReads);
    assert(basicData.qStart_.size()     == numReads);
    assert(basicData.qEnd_.size()       == numReads);
    assert(basicData.holeNumber_.size() == numReads);
    assert(basicData.readQual_.size()   == numReads);
    assert(basicData.ctxtFlag_.size()   == numReads);
    assert(basicData.fileOffset_.size() == numReads);
}

void PbiIndexIO::Save(const PbiRawData& index,
                      const std::string& filename)
{
    std::unique_ptr<BGZF, HtslibBgzfDeleter> bgzf(bgzf_open(filename.c_str(), "wb"));
    BGZF* fp = bgzf.get();
    if (fp == 0)
        throw std::runtime_error("could not open PBI file for writing");

    WriteHeader(index, fp);
    const uint32_t numReads = index.NumReads();
    if (numReads > 0) {
        WriteBasicData(index.BasicData(), numReads, fp);

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
    assert(barcodeData.bcForward_.size()   == numReads);
    assert(barcodeData.bcReverse_.size()   == numReads);
    assert(barcodeData.bcQual_.size()      == numReads);
    (void)numReads; // quash warnings building in release mode

    WriteBgzfVector(fp, barcodeData.bcForward_);
    WriteBgzfVector(fp, barcodeData.bcReverse_);
    WriteBgzfVector(fp, barcodeData.bcQual_);
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
    (void)numReads; // quash warnings building in release mode

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

void PbiIndexIO::WriteBasicData(const PbiRawBasicData& basicData,
                                const uint32_t numReads,
                                BGZF* fp)
{
    assert(basicData.rgId_.size()       == numReads);
    assert(basicData.qStart_.size()     == numReads);
    assert(basicData.qEnd_.size()       == numReads);
    assert(basicData.holeNumber_.size() == numReads);
    assert(basicData.readQual_.size()   == numReads);
    assert(basicData.ctxtFlag_.size()   == numReads);
    assert(basicData.fileOffset_.size() == numReads);
    (void)numReads; // quash warnings building in release mode

    WriteBgzfVector(fp, basicData.rgId_);
    WriteBgzfVector(fp, basicData.qStart_);
    WriteBgzfVector(fp, basicData.qEnd_);
    WriteBgzfVector(fp, basicData.holeNumber_);
    WriteBgzfVector(fp, basicData.readQual_);
    WriteBgzfVector(fp, basicData.ctxtFlag_);
    WriteBgzfVector(fp, basicData.fileOffset_);
}
