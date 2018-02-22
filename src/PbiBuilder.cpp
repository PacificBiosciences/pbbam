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

#include "PbbamInternalConfig.h"

#include "pbbam/PbiBuilder.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <thread>
#include <tuple>

#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>

#include <boost/numeric/conversion/cast.hpp>

#include "MemoryUtils.h"
#include "pbbam/BamRecord.h"
#include "pbbam/MakeUnique.h"
#include "pbbam/PbiRawData.h"

namespace {

using RecordType = PacBio::BAM::RecordType;
std::string ToString(const RecordType type)
{
    // clang-format off
    static const auto lookup = std::map<RecordType, std::string>
    {
        { RecordType::ZMW,        "ZMW" },
        { RecordType::HQREGION,   "HQREGION" },
        { RecordType::SUBREAD,    "SUBREAD" },
        { RecordType::CCS,        "CCS" },
        { RecordType::SCRAP,      "SCRAP" },
        { RecordType::UNKNOWN,    "UNKNOWN" }
    };
    // clang-format on

    try {
        return lookup.at(type);
    } catch (std::exception&) {
        throw std::runtime_error("error: unknown RecordType encountered");
    }
}

}  // namespace anonymous

namespace PacBio {
namespace BAM {
namespace internal {

template <typename T>
inline void SwapEndianness(std::vector<T>& data)
{
    const size_t elementSize = sizeof(T);
    const size_t numReads = data.size();
    switch (elementSize) {
        case 1:
            break;  // no swapping necessary
        case 2:
            for (size_t i = 0; i < numReads; ++i)
                ed_swap_2p(&data[i]);
            break;
        case 4:
            for (size_t i = 0; i < numReads; ++i)
                ed_swap_4p(&data[i]);
            break;
        case 8:
            for (size_t i = 0; i < numReads; ++i)
                ed_swap_8p(&data[i]);
            break;
        default:
            throw std::runtime_error("unsupported element size");
    }
}

void bgzf_write_safe(BGZF* fp, const void* data, size_t length)
{
    auto ret = ::bgzf_write(fp, data, length);
    if (ret < 0U) {
        auto msg = "Non-zero returned from bgzf_write(). Out of disk space?";
        throw std::runtime_error(msg);
    }
}

template <typename T>
inline void WriteBgzfVector(BGZF* fp, std::vector<T>& data, const size_t numElements)
{
    assert(fp);
    if (fp->is_be) SwapEndianness(data);
    bgzf_write_safe(fp, &data[0], numElements * sizeof(T));
}

// --------------------------
// PbiTempFile
// --------------------------

template <typename T>
class PbiTempFile
{
public:
    constexpr static size_t MaxBufferSize = 0x10000;  // 64K
    constexpr static size_t ElementSize = sizeof(T);
    constexpr static size_t MaxElementCount = MaxBufferSize / ElementSize;

public:
    PbiTempFile(std::string fn);
    ~PbiTempFile();

public:
    void Close();
    const std::vector<T>& Data() const;
    std::vector<T>& Data();
    void Flush();
    size_t Read(const size_t count);
    void Rewind();
    void Write(T value);

private:
    void WriteToFile();

private:
    // file info
    std::string fn_;
    std::unique_ptr<FILE, internal::FileDeleter> fp_ = nullptr;

    // data storage/tracking
    std::vector<T> buffer_;
    size_t numElementsWritten_ = 0;
};

template <typename T>
PbiTempFile<T>::PbiTempFile(std::string fn)
    : fn_{std::move(fn)}, fp_{std::fopen(fn_.c_str(), "w+b")}
{
    if (fp_ == nullptr) throw std::runtime_error("could not open temp file: " + fn_);
    buffer_.reserve(MaxElementCount);
}

template <typename T>
PbiTempFile<T>::~PbiTempFile()
{
    remove(fn_.c_str());
}

template <typename T>
void PbiTempFile<T>::Close()
{
    Flush();  // dtor will take care of closing file handle
}

template <typename T>
const std::vector<T>& PbiTempFile<T>::Data() const
{
    return buffer_;
}

template <typename T>
std::vector<T>& PbiTempFile<T>::Data()
{
    return buffer_;
}

template <typename T>
void PbiTempFile<T>::Flush()
{
    WriteToFile();
    buffer_.clear();
}

template <typename T>
size_t PbiTempFile<T>::Read(const size_t count)
{
    const auto actualCount = std::min(count, numElementsWritten_);
    buffer_.resize(actualCount);
    return fread(buffer_.data(), ElementSize, actualCount, fp_.get());
}

template <typename T>
void PbiTempFile<T>::Rewind()
{
    Flush();

    const auto ret = fseek(fp_.get(), 0, SEEK_SET);
    if (ret != 0) throw std::runtime_error("could not rewind temp file" + fn_);
}

template <typename T>
void PbiTempFile<T>::Write(T value)
{
    buffer_.push_back(value);

    // maybe flush
    if (buffer_.size() == MaxElementCount) Flush();
}

template <typename T>
void PbiTempFile<T>::WriteToFile()
{
    numElementsWritten_ += fwrite(buffer_.data(), ElementSize, buffer_.size(), fp_.get());
}

// --------------------------
// PbiReferenceDataBuilder
// --------------------------

class PbiReferenceDataBuilder
{
public:
    using ReferenceRows = std::pair<int32_t, int32_t>;  // [startRow, endRow)

public:
    explicit PbiReferenceDataBuilder(const size_t numReferenceSequences);

public:
    bool AddRecord(const BamRecord& record, const int32_t rowNumber);

    PbiRawReferenceData Result() const;

    void WriteData(BGZF* bgzf);

private:
    int32_t lastRefId_ = -1;
    Position lastPos_ = -1;
    std::map<uint32_t, PbiReferenceEntry> rawReferenceEntries_;
};

PbiReferenceDataBuilder::PbiReferenceDataBuilder(const size_t numReferenceSequences)
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

bool PbiReferenceDataBuilder::AddRecord(const BamRecord& record, const int32_t rowNumber)
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
                rawReferenceEntries_.at(PbiReferenceEntry::UNMAPPED_ID);
            if (unmappedEntry.beginRow_ != PbiReferenceEntry::UNSET_ROW) return false;

            // if we've already seen data for this new tId
            // (remember we're coming from another tId)
            //
            // error: refs are out of order (can stop checking refs)
            //
            PbiReferenceEntry& currentEntry = rawReferenceEntries_.at((uint32_t)tId);
            if (currentEntry.beginRow_ != PbiReferenceEntry::UNSET_ROW) return false;
        }
        lastRefId_ = tId;
    } else if (tId >= 0 && lastPos_ > pos)
        return false;  // error: positions out of order

    // update row numbers
    PbiReferenceEntry& entry = rawReferenceEntries_.at((uint32_t)tId);
    if (entry.beginRow_ == PbiReferenceEntry::UNSET_ROW) entry.beginRow_ = rowNumber;
    entry.endRow_ = rowNumber + 1;

    // update pos (for sorting check next go-round)
    lastPos_ = pos;
    return true;
}

PbiRawReferenceData PbiReferenceDataBuilder::Result() const
{
    // PbiReferenceEntries will be sorted thanks to std::map
    // tId will be at end since we're sorting on the uint cast of -1
    PbiRawReferenceData result;
    result.entries_.reserve(rawReferenceEntries_.size());
    for (auto& entry : rawReferenceEntries_)
        result.entries_.push_back(entry.second);
    return result;
}

void PbiReferenceDataBuilder::WriteData(BGZF* bgzf)
{
    const auto refData = Result();

    // num_refs
    uint32_t numRefs = refData.entries_.size();
    if (bgzf->is_be) numRefs = ed_swap_4(numRefs);
    bgzf_write_safe(bgzf, &numRefs, 4);

    // reference entries
    numRefs = refData.entries_.size();  // need to reset after maybe endian-swapping
    for (size_t i = 0; i < numRefs; ++i) {
        auto& entry = refData.entries_[i];
        auto tId = entry.tId_;
        auto beginRow = entry.beginRow_;
        auto endRow = entry.endRow_;
        if (bgzf->is_be) {
            tId = ed_swap_4(tId);
            beginRow = ed_swap_4(beginRow);
            endRow = ed_swap_4(endRow);
        }
        bgzf_write_safe(bgzf, &tId, 4);
        bgzf_write_safe(bgzf, &beginRow, 4);
        bgzf_write_safe(bgzf, &endRow, 4);
    }
}

// --------------------------------------------
// PbiBuilderPrivate - builder implementation
// --------------------------------------------

class PbiBuilderPrivate
{
    using CompressionLevel = PacBio::BAM::PbiBuilder::CompressionLevel;

public:
    PbiBuilderPrivate(const std::string& pbiFilename, const size_t numReferenceSequences,
                      const bool isCoordinateSorted, const CompressionLevel compressionLevel,
                      const size_t numThreads);

    ~PbiBuilderPrivate() noexcept;

public:
    void AddRecord(const BamRecord& record, const int64_t vOffset);
    void Close();

private:
    // store record data
    void AddBarcodeData(const BamRecord& record);
    void AddBasicData(const BamRecord& record, const int64_t vOffset);
    void AddMappedData(const BamRecord& record);
    void AddReferenceData(const BamRecord& record, const uint32_t currentRow);

    // read from temp files & write PBI data
    void OpenPbiFile();

    template <typename T>
    void WriteFromTempFile(PbiTempFile<T>& tempFile, BGZF* bgzf);

    void WritePbiHeader(BGZF* bgzf);
    void WriteReferenceData(BGZF* bgzf);

private:
    // basic data
    PbiTempFile<int32_t> rgIdFile_;
    PbiTempFile<int32_t> qStartFile_;
    PbiTempFile<int32_t> qEndFile_;
    PbiTempFile<int32_t> holeNumFile_;
    PbiTempFile<float> readQualFile_;
    PbiTempFile<uint8_t> ctxtFile_;
    PbiTempFile<int64_t> fileOffsetFile_;

    // mapped data
    PbiTempFile<int32_t> tIdFile_;
    PbiTempFile<uint32_t> tStartFile_;
    PbiTempFile<uint32_t> tEndFile_;
    PbiTempFile<uint32_t> aStartFile_;
    PbiTempFile<uint32_t> aEndFile_;
    PbiTempFile<uint8_t> revStrandFile_;
    PbiTempFile<uint32_t> nMFile_;
    PbiTempFile<uint32_t> nMMFile_;
    PbiTempFile<uint8_t> mapQualFile_;

    // barcode data
    PbiTempFile<int16_t> bcForwardFile_;
    PbiTempFile<int16_t> bcReverseFile_;
    PbiTempFile<int8_t> bcQualFile_;

    // reference data
    std::unique_ptr<PbiReferenceDataBuilder> refDataBuilder_ = nullptr;

    // output file info
    std::string pbiFilename_;
    std::unique_ptr<BGZF, internal::HtslibBgzfDeleter> outFile_ = nullptr;
    CompressionLevel compressionLevel_;
    size_t numThreads_;

    // tracking data
    uint32_t currentRow_ = 0;
    bool isClosed_ = false;
    bool hasBarcodeData_ = false;
    bool hasMappedData_ = false;
};

PbiBuilderPrivate::PbiBuilderPrivate(const std::string& pbiFilename,
                                     const size_t numReferenceSequences,
                                     const bool isCoordinateSorted,
                                     const CompressionLevel compressionLevel,
                                     const size_t numThreads)
    : rgIdFile_{pbiFilename + ".rgId.tmp"}
    , qStartFile_{pbiFilename + ".qStart.tmp"}
    , qEndFile_{pbiFilename + ".qEnd.tmp"}
    , holeNumFile_{pbiFilename + ".holeNum.tmp"}
    , readQualFile_{pbiFilename + ".rq.tmp"}
    , ctxtFile_{pbiFilename + ".ctxt.tmp"}
    , fileOffsetFile_{pbiFilename + ".offset.tmp"}
    , tIdFile_{pbiFilename + ".tId.tmp"}
    , tStartFile_{pbiFilename + ".tStart.tmp"}
    , tEndFile_{pbiFilename + ".tEnd.tmp"}
    , aStartFile_{pbiFilename + ".aStart.tmp"}
    , aEndFile_{pbiFilename + ".aEnd.tmp"}
    , revStrandFile_{pbiFilename + ".revStrand.tmp"}
    , nMFile_{pbiFilename + ".nm.tmp"}
    , nMMFile_{pbiFilename + ".nmm.tmp"}
    , mapQualFile_{pbiFilename + ".mapQual.tmp"}
    , bcForwardFile_{pbiFilename + ".bcForward.tmp"}
    , bcReverseFile_{pbiFilename + ".bcReverse.tmp"}
    , bcQualFile_{pbiFilename + ".bcQual.tmp"}
    , pbiFilename_{pbiFilename}
    , compressionLevel_{compressionLevel}
    , numThreads_{numThreads}
{
    if (isCoordinateSorted && numReferenceSequences > 0)
        refDataBuilder_ = std::make_unique<PbiReferenceDataBuilder>(numReferenceSequences);
}

PbiBuilderPrivate::~PbiBuilderPrivate() noexcept
{
    if (!isClosed_) {
        try {
            Close();
        } catch (...) {
            // swallow any exceptions & remain no-throw from dtor
        }
    }
}

void PbiBuilderPrivate::AddBarcodeData(const BamRecord& b)
{
    // initialize w/ 'missing' value
    auto bcForward = int16_t{-1};
    auto bcReverse = int16_t{-1};
    auto bcQuality = int8_t{-1};

    // check for any barcode data (both required)
    if (b.HasBarcodes() && b.HasBarcodeQuality()) {
        // fetch data from record
        std::tie(bcForward, bcReverse) = b.Barcodes();
        bcQuality = static_cast<int8_t>(b.BarcodeQuality());

        // double-check & reset to 'missing' value if any less than zero
        if (bcForward < 0 && bcReverse < 0 && bcQuality < 0) {
            bcForward = -1;
            bcReverse = -1;
            bcQuality = -1;
        } else
            hasBarcodeData_ = true;
    }

    // store
    bcForwardFile_.Write(bcForward);
    bcReverseFile_.Write(bcReverse);
    bcQualFile_.Write(bcQuality);
}

void PbiBuilderPrivate::AddBasicData(const BamRecord& b, const int64_t vOffset)
{
    // read group ID
    const auto rgId = [&b]() -> int32_t {
        auto rgIdString = b.ReadGroupId();
        if (rgIdString.empty()) rgIdString = MakeReadGroupId(b.MovieName(), ToString(b.Type()));
        const auto rawId = std::stoul(rgIdString, nullptr, 16);
        return static_cast<int32_t>(rawId);
    }();

    // query start/end
    const auto isCcs = (b.Type() == RecordType::CCS);
    const auto qStart = int32_t{(isCcs ? -1 : b.QueryStart())};
    const auto qEnd = int32_t{(isCcs ? -1 : b.QueryEnd())};

    // add'l data
    const auto holeNum = int32_t{(b.HasHoleNumber() ? b.HoleNumber() : 0)};
    const auto readAccuracy =
        float{(b.HasReadAccuracy() ? boost::numeric_cast<float>(b.ReadAccuracy()) : 0.0F)};
    const auto ctxt = uint8_t{
        (b.HasLocalContextFlags() ? b.LocalContextFlags() : LocalContextFlags::NO_LOCAL_CONTEXT)};

    // store
    rgIdFile_.Write(rgId);
    qStartFile_.Write(qStart);
    qEndFile_.Write(qEnd);
    holeNumFile_.Write(holeNum);
    ctxtFile_.Write(ctxt);
    readQualFile_.Write(readAccuracy);
    fileOffsetFile_.Write(vOffset);
}

void PbiBuilderPrivate::AddMappedData(const BamRecord& b)
{
    // fetch data
    const auto tId = int32_t{b.ReferenceId()};
    const auto tStart = static_cast<uint32_t>(b.ReferenceStart());
    const auto tEnd = static_cast<uint32_t>(b.ReferenceEnd());
    const auto aStart = static_cast<uint32_t>(b.AlignedStart());
    const auto aEnd = static_cast<uint32_t>(b.AlignedEnd());

    const auto isReverseStrand = [&b]() -> uint8_t {
        return (b.AlignedStrand() == Strand::REVERSE ? 1 : 0);
    }();

    const auto matchData = b.NumMatchesAndMismatches();
    const auto nM = static_cast<uint32_t>(matchData.first);
    const auto nMM = static_cast<uint32_t>(matchData.second);
    const auto mapQuality = uint8_t{b.MapQuality()};

    if (tId >= 0) hasMappedData_ = true;

    // store
    tIdFile_.Write(tId);
    tStartFile_.Write(tStart);
    tEndFile_.Write(tEnd);
    aStartFile_.Write(aStart);
    aEndFile_.Write(aEnd);
    revStrandFile_.Write(isReverseStrand);
    nMFile_.Write(nM);
    nMMFile_.Write(nMM);
    mapQualFile_.Write(mapQuality);
}

void PbiBuilderPrivate::AddRecord(const BamRecord& b, const int64_t vOffset)
{
    // ensure updated data
    internal::BamRecordMemory::UpdateRecordTags(b);
    b.ResetCachedPositions();

    // store data
    AddBasicData(b, vOffset);
    AddMappedData(b);
    AddBarcodeData(b);
    AddReferenceData(b, currentRow_);

    // increment row counter
    ++currentRow_;
}

void PbiBuilderPrivate::AddReferenceData(const BamRecord& b, const uint32_t currentRow)
{
    // only add if coordinate-sorted hint is set
    // update with info from refDataBuilder
    if (refDataBuilder_) {
        const auto sorted = refDataBuilder_->AddRecord(b, currentRow);
        if (!sorted) refDataBuilder_.reset(nullptr);
    }
}

void PbiBuilderPrivate::Close()
{
    // open PBI file for writing
    OpenPbiFile();
    auto* bgzf = outFile_.get();

    // header section
    WritePbiHeader(bgzf);

    // 'basic' data section
    WriteFromTempFile(rgIdFile_, bgzf);
    WriteFromTempFile(qStartFile_, bgzf);
    WriteFromTempFile(qEndFile_, bgzf);
    WriteFromTempFile(holeNumFile_, bgzf);
    WriteFromTempFile(readQualFile_, bgzf);
    WriteFromTempFile(ctxtFile_, bgzf);
    WriteFromTempFile(fileOffsetFile_, bgzf);

    // mapped data section
    if (hasMappedData_) {
        WriteFromTempFile(tIdFile_, bgzf);
        WriteFromTempFile(tStartFile_, bgzf);
        WriteFromTempFile(tEndFile_, bgzf);
        WriteFromTempFile(aStartFile_, bgzf);
        WriteFromTempFile(aEndFile_, bgzf);
        WriteFromTempFile(revStrandFile_, bgzf);
        WriteFromTempFile(nMFile_, bgzf);
        WriteFromTempFile(nMMFile_, bgzf);
        WriteFromTempFile(mapQualFile_, bgzf);
    }

    // reference data section
    if (refDataBuilder_) WriteReferenceData(bgzf);

    // barcode data section
    if (hasBarcodeData_) {
        WriteFromTempFile(bcForwardFile_, bgzf);
        WriteFromTempFile(bcReverseFile_, bgzf);
        WriteFromTempFile(bcQualFile_, bgzf);
    }

    // finally, set flag
    isClosed_ = true;
}

void PbiBuilderPrivate::OpenPbiFile()
{
    // open file handle
    const auto mode = std::string("wb") + std::to_string(static_cast<int>(compressionLevel_));
    outFile_.reset(bgzf_open(pbiFilename_.c_str(), mode.c_str()));
    if (outFile_ == nullptr) throw std::runtime_error("could not open output file");

    // if no explicit thread count given, attempt built-in check
    size_t actualNumThreads = numThreads_;
    if (actualNumThreads == 0) {
        actualNumThreads = std::thread::hardware_concurrency();

        // if still unknown, default to single-threaded
        if (actualNumThreads == 0) actualNumThreads = 1;
    }

    // if multithreading requested, enable it
    if (actualNumThreads > 1) bgzf_mt(outFile_.get(), actualNumThreads, 256);
}

template <typename T>
void PbiBuilderPrivate::WriteFromTempFile(PbiTempFile<T>& tempFile, BGZF* bgzf)
{
    using TempFileType = PbiTempFile<T>;
    static constexpr const auto maxElementCount = TempFileType::MaxElementCount;

    tempFile.Rewind();

    size_t totalNumRead = 0;
    for (size_t i = 0; totalNumRead < currentRow_; ++i) {
        const auto numRead = tempFile.Read(maxElementCount);
        auto& data = tempFile.Data();
        WriteBgzfVector(bgzf, data, numRead);
        totalNumRead += numRead;
    }
}

void PbiBuilderPrivate::WritePbiHeader(BGZF* bgzf)
{
    // 'magic' string
    static constexpr const std::array<char, 4> magic{{'P', 'B', 'I', '\1'}};
    bgzf_write_safe(bgzf, magic.data(), 4);

    PbiFile::Sections sections = PbiFile::BASIC;
    if (hasMappedData_) sections |= PbiFile::MAPPED;
    if (hasBarcodeData_) sections |= PbiFile::BARCODE;
    if (refDataBuilder_) sections |= PbiFile::REFERENCE;

    // version, pbi_flags, & n_reads
    auto version = static_cast<uint32_t>(PbiFile::CurrentVersion);
    auto pbi_flags = static_cast<uint16_t>(sections);
    auto numReads = currentRow_;
    if (bgzf->is_be) {
        version = ed_swap_4(version);
        pbi_flags = ed_swap_2(pbi_flags);
        numReads = ed_swap_4(numReads);
    }
    bgzf_write_safe(bgzf, &version, 4);
    bgzf_write_safe(bgzf, &pbi_flags, 2);
    bgzf_write_safe(bgzf, &numReads, 4);

    // reserved space
    char reserved[18];
    memset(reserved, 0, 18);
    bgzf_write_safe(bgzf, reserved, 18);
}

void PbiBuilderPrivate::WriteReferenceData(BGZF* bgzf) { refDataBuilder_->WriteData(bgzf); }

}  // namespace internal

// --------------------------------------------
// PbiBuilder - builder API
// --------------------------------------------

PbiBuilder::PbiBuilder(const std::string& pbiFilename, const CompressionLevel compressionLevel,
                       const size_t numThreads)
    : PbiBuilder{pbiFilename, 0, false, compressionLevel, numThreads}
{
}

PbiBuilder::PbiBuilder(const std::string& pbiFilename, const size_t numReferenceSequences,
                       const CompressionLevel compressionLevel, const size_t numThreads)
    : PbiBuilder{pbiFilename, numReferenceSequences, (numReferenceSequences > 0), compressionLevel,
                 numThreads}
{
}

PbiBuilder::PbiBuilder(const std::string& pbiFilename, const size_t numReferenceSequences,
                       const bool isCoordinateSorted, const CompressionLevel compressionLevel,
                       const size_t numThreads)
    : d_{new internal::PbiBuilderPrivate{pbiFilename, numReferenceSequences, isCoordinateSorted,
                                         compressionLevel, numThreads}}
{
}

PbiBuilder::~PbiBuilder() noexcept {}

void PbiBuilder::AddRecord(const BamRecord& record, const int64_t vOffset)
{
    d_->AddRecord(record, vOffset);
}

void PbiBuilder::Close() { d_->Close(); }

}  // namespace BAM
}  // namespace PacBio
