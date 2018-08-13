/// File Description
/// \file IndexedBamWriter.cpp
/// \brief Implements the IndexedBamWriter class
//
// Author: Derek Barnett

#include "pbbam/IndexedBamWriter.h"

#include <sys/stat.h>

#include <array>
#include <atomic>
#include <cassert>
#include <condition_variable>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <thread>

#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>

#include "pbbam/BamHeader.h"
#include "pbbam/BamRecord.h"
#include "pbbam/BamWriter.h"
#include "pbbam/MakeUnique.h"
#include "pbbam/PbiRawData.h"
#include "pbbam/Unused.h"
#include "pbbam/Validator.h"

#include "FileProducer.h"
#include "MemoryUtils.h"

namespace PacBio {
namespace BAM {
namespace internal {

void bgzf_write_safe2(BGZF* fp, const void* data, size_t length)
{
    const auto ret = bgzf_write(fp, data, length);
    if (ret < 0L)
        throw std::runtime_error{"Non-zero returned from bgzf_write(). Out of disk space?"};
}

struct GzIndexEntry
{
    int64_t vAddress;
    int64_t uAddress;
};

template <typename T>
inline void SwapEndianness2(std::vector<T>& data)
{
    constexpr const size_t elementSize = sizeof(T);
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
            throw std::runtime_error{"unsupported element size"};
    }
}

template <typename T>
inline void WriteBgzfVector2(BGZF* fp, std::vector<T>& data, const size_t numElements)
{
    assert(fp);
    if (fp->is_be) SwapEndianness2(data);
    bgzf_write_safe2(fp, data.data(), numElements * sizeof(T));
}

using RecordType = PacBio::BAM::RecordType;
std::string ToString2(const RecordType type)
{
    // clang-format off
    static const auto lookup = std::map<RecordType, std::string>
    {
        { RecordType::ZMW,        "ZMW" },
        { RecordType::HQREGION,   "HQREGION" },
        { RecordType::SUBREAD,    "SUBREAD" },
        { RecordType::CCS,        "CCS" },
        { RecordType::SCRAP,      "SCRAP" },
        { RecordType::TRANSCRIPT, "TRANSCRIPT" },
        { RecordType::UNKNOWN,    "UNKNOWN" }
    };
    // clang-format on

    try {
        return lookup.at(type);
    } catch (const std::exception&) {
        throw std::runtime_error{"error: unknown RecordType encountered"};
    }
}

template <typename T>
class PbiTempFile2
{
public:
    constexpr static const size_t ElementSize = sizeof(T);

public:
    PbiTempFile2(std::string fn, size_t maxBufferSize)
        : fn_{std::move(fn)}
        , fp_{std::fopen(fn_.c_str(), "w+b")}
        , maxElementCount_{maxBufferSize / ElementSize}
    {
        if (fp_ == nullptr) throw std::runtime_error{"could not open temp file: " + fn_};
        buffer_.reserve(maxElementCount_);
    }
    ~PbiTempFile2() { remove(fn_.c_str()); }

public:
    void Close()
    {
        Flush();  // dtor will take care of closing file handle
    }

    const std::vector<T>& Data() const { return buffer_; }

    std::vector<T>& Data() { return buffer_; }

    void Flush()
    {
        WriteToFile();
        buffer_.clear();
    }

    size_t Read(const size_t count)
    {
        const auto actualCount = std::min(count, numElementsWritten_);
        buffer_.resize(actualCount);
        return fread(buffer_.data(), ElementSize, actualCount, fp_.get());
    }

    void Rewind()
    {
        Flush();

        const auto ret = fseek(fp_.get(), 0, SEEK_SET);
        if (ret != 0) throw std::runtime_error{"could not rewind temp file" + fn_};
    }

    void Write(T value)
    {
        buffer_.push_back(value);

        // maybe flush
        if (buffer_.size() == maxElementCount_) Flush();
    }

    size_t MaxElementCount() const { return maxElementCount_; }

private:
    void WriteToFile()
    {
        numElementsWritten_ += fwrite(buffer_.data(), ElementSize, buffer_.size(), fp_.get());
    }

private:
    // file info
    std::string fn_;
    std::unique_ptr<FILE, internal::FileDeleter> fp_;

    // data storage/tracking
    std::vector<T> buffer_;
    size_t numElementsWritten_ = 0;
    size_t maxElementCount_;
};

class PbiReferenceDataBuilder2
{
public:
    using ReferenceRows = std::pair<int32_t, int32_t>;  // [startRow, endRow)

public:
    explicit PbiReferenceDataBuilder2(const size_t numReferenceSequences)
    {
        // initialize with number of references we expect to see
        //
        // we can add more later, but want to ensure known references have an entry
        // even if no records are observed mapping to it
        //
        for (size_t i = 0; i < numReferenceSequences; ++i)
            rawReferenceEntries_[i] = PbiReferenceEntry(i);

        // also create an "unmapped" entry
        rawReferenceEntries_[PbiReferenceEntry::UNMAPPED_ID] = PbiReferenceEntry{};
    }

public:
    bool AddRecord(const BamRecord& record, const int32_t rowNumber)
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
                PbiReferenceEntry& currentEntry =
                    rawReferenceEntries_.at(static_cast<uint32_t>(tId));
                if (currentEntry.beginRow_ != PbiReferenceEntry::UNSET_ROW) return false;
            }
            lastRefId_ = tId;
        } else if (tId >= 0 && lastPos_ > pos)
            return false;  // error: positions out of order

        // update row numbers
        PbiReferenceEntry& entry = rawReferenceEntries_.at(static_cast<uint32_t>(tId));
        if (entry.beginRow_ == PbiReferenceEntry::UNSET_ROW) entry.beginRow_ = rowNumber;
        entry.endRow_ = rowNumber + 1;

        // update pos (for sorting check next go-round)
        lastPos_ = pos;
        return true;
    }

    PbiRawReferenceData Result() const
    {
        // PbiReferenceEntries will be sorted thanks to std::map
        // tId will be at end since we're sorting on the uint cast of -1
        PbiRawReferenceData result;
        result.entries_.reserve(rawReferenceEntries_.size());
        for (const auto& entry : rawReferenceEntries_)
            result.entries_.push_back(entry.second);
        return result;
    }

    void WriteData(BGZF* bgzf)
    {
        const auto refData = Result();

        // num_refs
        uint32_t numRefs = refData.entries_.size();
        if (bgzf->is_be) numRefs = ed_swap_4(numRefs);
        bgzf_write_safe2(bgzf, &numRefs, 4);

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
            bgzf_write_safe2(bgzf, &tId, 4);
            bgzf_write_safe2(bgzf, &beginRow, 4);
            bgzf_write_safe2(bgzf, &endRow, 4);
        }
    }

private:
    int32_t lastRefId_ = -1;
    Position lastPos_ = -1;
    std::map<uint32_t, PbiReferenceEntry> rawReferenceEntries_;
};

// TODO: come back to refseqs, sorting, etc
class PbiBuilder2
{
public:
    PbiBuilder2(const std::string& bamFilename, const std::string& pbiFilename,
                const PbiBuilder::CompressionLevel compressionLevel, const size_t numThreads,
                const size_t fileBufferSize)
        //                const size_t numReferenceSequences = 0
        //                const bool isCoordinateSorted = false
        : rgIdFile_{pbiFilename + ".rgId.tmp", fileBufferSize},
          qStartFile_{pbiFilename + ".qStart.tmp", fileBufferSize},
          qEndFile_{pbiFilename + ".qEnd.tmp", fileBufferSize},
          holeNumFile_{pbiFilename + ".holeNum.tmp", fileBufferSize},
          readQualFile_{pbiFilename + ".rq.tmp", fileBufferSize},
          ctxtFile_{pbiFilename + ".ctxt.tmp", fileBufferSize},
          fileOffsetFile_{pbiFilename + ".offset.tmp", fileBufferSize},
          tIdFile_{pbiFilename + ".tId.tmp", fileBufferSize},
          tStartFile_{pbiFilename + ".tStart.tmp", fileBufferSize},
          tEndFile_{pbiFilename + ".tEnd.tmp", fileBufferSize},
          aStartFile_{pbiFilename + ".aStart.tmp", fileBufferSize},
          aEndFile_{pbiFilename + ".aEnd.tmp", fileBufferSize},
          revStrandFile_{pbiFilename + ".revStrand.tmp", fileBufferSize},
          nMFile_{pbiFilename + ".nm.tmp", fileBufferSize},
          nMMFile_{pbiFilename + ".nmm.tmp", fileBufferSize},
          mapQualFile_{pbiFilename + ".mapQual.tmp", fileBufferSize},
          bcForwardFile_{pbiFilename + ".bcForward.tmp", fileBufferSize},
          bcReverseFile_{pbiFilename + ".bcReverse.tmp", fileBufferSize},
          bcQualFile_{pbiFilename + ".bcQual.tmp", fileBufferSize},
          bamFilename_{bamFilename},
          pbiFilename_{pbiFilename},
          compressionLevel_{compressionLevel},
          numThreads_{numThreads}
    {
        // TODO: come back to this
        //        if (isCoordinateSorted && numReferenceSequences > 0)
        //            refDataBuilder_ = std::make_unique<PbiReferenceDataBuilder2>(numReferenceSequences);
    }

    void AddRecord(const BamRecord& b, const int64_t uOffset)
    {
        // ensure updated data (?)
        internal::BamRecordMemory::UpdateRecordTags(b);
        b.ResetCachedPositions();

        // store data
        AddBasicData(b, uOffset);
        AddMappedData(b);
        AddBarcodeData(b);
        AddReferenceData(b, currentRow_);

        // increment row counter
        ++currentRow_;
    }

    void Close()
    {
        if (isClosed_) return;

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

        // special vOffset transformation
        WriteVirtualOffsets(bgzf);

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

private:
    // store record data
    void AddBarcodeData(const BamRecord& b)
    {
        // initialize w/ 'missing' value
        int16_t bcForward = -1;
        int16_t bcReverse = -1;
        int8_t bcQuality = -1;

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

    void AddBasicData(const BamRecord& b, const int64_t vOffset)
    {
        // read group ID
        const auto rgId = [&b]() -> int32_t {
            auto rgIdString = b.ReadGroupId();
            if (rgIdString.empty())
                rgIdString = MakeReadGroupId(b.MovieName(), ToString2(b.Type()));
            return std::stoul(rgIdString, nullptr, 16);
        }();

        // query start/end
        const auto isCcsOrTranscript = (IsCcsOrTranscript(b.Type()));
        const int32_t qStart = (isCcsOrTranscript ? -1 : b.QueryStart());
        const int32_t qEnd = (isCcsOrTranscript ? -1 : b.QueryEnd());

        // add'l data
        const int32_t holeNum = (b.HasHoleNumber() ? b.HoleNumber() : 0);
        const float readAccuracy =
            (b.HasReadAccuracy() ? boost::numeric_cast<float>(b.ReadAccuracy()) : 0.0F);
        const uint8_t ctxt = (b.HasLocalContextFlags() ? b.LocalContextFlags()
                                                       : LocalContextFlags::NO_LOCAL_CONTEXT);

        // store
        rgIdFile_.Write(rgId);
        qStartFile_.Write(qStart);
        qEndFile_.Write(qEnd);
        holeNumFile_.Write(holeNum);
        ctxtFile_.Write(ctxt);
        readQualFile_.Write(readAccuracy);
        fileOffsetFile_.Write(vOffset);
    }

    void AddMappedData(const BamRecord& b)
    {
        // fetch data
        const auto tId = b.ReferenceId();
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
        const auto mapQuality = b.MapQuality();

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

    void AddReferenceData(const BamRecord& b, const uint32_t currentRow)
    {
        // only add if coordinate-sorted hint is set
        // update with info from refDataBuilder
        if (refDataBuilder_) {
            const auto sorted = refDataBuilder_->AddRecord(b, currentRow);
            if (!sorted) refDataBuilder_.reset(nullptr);
        }
    }

    // read from temp files & write PBI data
    void OpenPbiFile()
    {
        // open file handle
        const auto mode = std::string("wb") + std::to_string(static_cast<int>(compressionLevel_));
        outFile_.reset(bgzf_open(pbiFilename_.c_str(), mode.c_str()));
        if (outFile_ == nullptr) throw std::runtime_error{"could not open output file"};

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
    void WriteFromTempFile(PbiTempFile2<T>& tempFile, BGZF* bgzf)
    {
        const auto maxElementCount = tempFile.MaxElementCount();

        tempFile.Rewind();

        size_t totalNumRead = 0;
        for (size_t i = 0; totalNumRead < currentRow_; ++i) {
            const auto numRead = tempFile.Read(maxElementCount);
            auto& data = tempFile.Data();
            WriteBgzfVector2(bgzf, data, numRead);
            totalNumRead += numRead;
        }
    }

    void WritePbiHeader(BGZF* bgzf)
    {
        // 'magic' string
        static constexpr const std::array<char, 4> magic{{'P', 'B', 'I', '\1'}};
        bgzf_write_safe2(bgzf, magic.data(), 4);

        PbiFile::Sections sections = PbiFile::BASIC;
        if (hasMappedData_) sections |= PbiFile::MAPPED;
        if (hasBarcodeData_) sections |= PbiFile::BARCODE;
        if (refDataBuilder_) sections |= PbiFile::REFERENCE;

        // version, pbi_flags, & n_reads
        auto version = static_cast<uint32_t>(PbiFile::CurrentVersion);
        uint16_t pbi_flags = sections;
        auto numReads = currentRow_;
        if (bgzf->is_be) {
            version = ed_swap_4(version);
            pbi_flags = ed_swap_2(pbi_flags);
            numReads = ed_swap_4(numReads);
        }
        bgzf_write_safe2(bgzf, &version, 4);
        bgzf_write_safe2(bgzf, &pbi_flags, 2);
        bgzf_write_safe2(bgzf, &numReads, 4);

        // reserved space
        char reserved[18];
        memset(reserved, 0, 18);
        bgzf_write_safe2(bgzf, reserved, 18);
    }

    void WriteReferenceData(BGZF* bgzf) { refDataBuilder_->WriteData(bgzf); }

    std::vector<GzIndexEntry> LoadGzi()
    {
        //
        // Open GZI file & load its contents. About to use for offset transformation.
        //

        const std::string gziFn{bamFilename_ + ".gzi"};
        std::unique_ptr<FILE, internal::FileDeleter> gziFile{fopen(gziFn.c_str(), "rb")};
        if (!gziFile) throw std::runtime_error{"could not open gzi file"};

        uint64_t numElements;
        const auto ret = fread(&numElements, sizeof(numElements), 1, gziFile.get());
        if (ret != 1) throw std::runtime_error{"could not read from gziFile"};
        if (ed_is_big()) ed_swap_8(numElements);

        std::vector<GzIndexEntry> result;
        result.reserve(numElements);
        for (uint32_t i = 0; i < numElements; ++i) {
            GzIndexEntry entry;
            fread(&entry.vAddress, sizeof(entry.vAddress), 1, gziFile.get());
            fread(&entry.uAddress, sizeof(entry.uAddress), 1, gziFile.get());
            if (ed_is_big()) {
                ed_swap_8(entry.vAddress);
                ed_swap_8(entry.uAddress);
            }
            result.push_back(std::move(entry));
        }
        return result;
    }

    void WriteVirtualOffsets(BGZF* bgzf)
    {
        auto index = LoadGzi();
        if (index.empty()) throw std::runtime_error{"empty GZI file"};
        std::sort(index.begin(), index.end(),
                  [](const GzIndexEntry& lhs, const GzIndexEntry& rhs) -> bool {
                      return lhs.uAddress < rhs.uAddress;
                  });

        const auto maxElementCount = fileOffsetFile_.MaxElementCount();
        fileOffsetFile_.Rewind();

        size_t k = 0;
        size_t totalNumRead = 0;
        for (size_t i = 0; totalNumRead < currentRow_; ++i) {

            const auto numRead = fileOffsetFile_.Read(maxElementCount);
            auto& data = fileOffsetFile_.Data();
            for (size_t j = 0; j < data.size(); ++j) {
                if (k < index.size() - 1 &&
                    (static_cast<uint64_t>(index.at(k + 1).uAddress) <= data[j]))
                    ++k;
                const GzIndexEntry& e = index.at(k);
                const int64_t uOffset = data[j] - e.uAddress;
                const auto result = ((e.vAddress << 16) | uOffset);
                data[j] = result;
            }
            WriteBgzfVector2(bgzf, data, numRead);
            totalNumRead += numRead;
        }
    }

private:
    // basic data
    PbiTempFile2<int32_t> rgIdFile_;
    PbiTempFile2<int32_t> qStartFile_;
    PbiTempFile2<int32_t> qEndFile_;
    PbiTempFile2<int32_t> holeNumFile_;
    PbiTempFile2<float> readQualFile_;
    PbiTempFile2<uint8_t> ctxtFile_;
    PbiTempFile2<uint64_t> fileOffsetFile_;

    // mapped data
    PbiTempFile2<int32_t> tIdFile_;
    PbiTempFile2<uint32_t> tStartFile_;
    PbiTempFile2<uint32_t> tEndFile_;
    PbiTempFile2<uint32_t> aStartFile_;
    PbiTempFile2<uint32_t> aEndFile_;
    PbiTempFile2<uint8_t> revStrandFile_;
    PbiTempFile2<uint32_t> nMFile_;
    PbiTempFile2<uint32_t> nMMFile_;
    PbiTempFile2<uint8_t> mapQualFile_;

    // barcode data
    PbiTempFile2<int16_t> bcForwardFile_;
    PbiTempFile2<int16_t> bcReverseFile_;
    PbiTempFile2<int8_t> bcQualFile_;

    // reference data
    std::unique_ptr<PbiReferenceDataBuilder2> refDataBuilder_ = nullptr;

    // output file info
    std::string bamFilename_;
    std::string pbiFilename_;
    std::unique_ptr<BGZF, internal::HtslibBgzfDeleter> outFile_ = nullptr;
    PbiBuilder::CompressionLevel compressionLevel_;
    size_t numThreads_;

    // tracking data
    uint32_t currentRow_ = 0;
    bool isClosed_ = false;
    bool hasBarcodeData_ = false;
    bool hasMappedData_ = false;
};

class IndexedBamWriterPrivate2  //: public internal::FileProducer
{

public:
    IndexedBamWriterPrivate2(const std::string& outputFilename, std::shared_ptr<bam_hdr_t> header,
                             const BamWriter::CompressionLevel bamCompressionLevel,
                             const size_t numBamThreads,
                             const PbiBuilder::CompressionLevel pbiCompressionLevel,
                             const size_t numPbiThreads, const size_t numGziThreads,
                             const size_t tempFileBufferSize)
        : bamFilename_{outputFilename}, header_{header}
    {
        OpenBam(bamCompressionLevel, numBamThreads);
        OpenGzi(numGziThreads);
        OpenPbi(pbiCompressionLevel, numPbiThreads, tempFileBufferSize);

        isOpen_ = true;
    }

    ~IndexedBamWriterPrivate2() noexcept
    {
        if (isOpen_) {
            try {
                Close();
            } catch (...) {
                // swallow any exceptions & remain no-throw from dtor
            }
        }
    }

    void Close()
    {
        // NOTE: keep this order of closing ( BAM -> GZI -> PBI )
        CloseBam();
        CloseGzi();
        ClosePbi();

        isOpen_ = false;
    }

    void CloseBam()
    {
        const auto ret = bgzf_flush(bam_.get()->fp.bgzf);
        UNUSED(ret);
        bam_.reset();
    }

    void CloseGzi()
    {
        done_ = true;
        gziThread_.join();
    }

    void ClosePbi() { builder_->Close(); }

    void OpenBam(const BamWriter::CompressionLevel compressionLevel, const size_t numThreads)
    {
        //
        // TODO: Compression level & numThreads are hardcoded here. Ok for
        //       prototyping but need to be tune-able via API.
        //

        if (!header_) throw std::runtime_error{"null header"};

        // open output BAM
        const auto usingFilename = bamFilename_;  //= TempFilename();
        const auto mode = std::string("wb") + std::to_string(static_cast<int>(compressionLevel));
        bam_.reset(sam_open(usingFilename.c_str(), mode.c_str()));
        if (!bam_) throw std::runtime_error{"could not open file for writing"};

        // maybe set multithreaded writing
        size_t actualNumThreads = numThreads;
        if (actualNumThreads == 0) {
            actualNumThreads = std::thread::hardware_concurrency();

            // if still unknown, default to single-threaded
            if (actualNumThreads == 0) actualNumThreads = 1;
        }
        if (actualNumThreads > 1) hts_set_threads(bam_.get(), actualNumThreads);

        // write header
        auto ret = sam_hdr_write(bam_.get(), header_.get());
        if (ret != 0) throw std::runtime_error{"could not write header"};
        ret = bgzf_flush(bam_.get()->fp.bgzf);

        // store file positions after header
        auto headerLength = [](const bam_hdr_t* hdr) -> size_t {
            const size_t textHeader = 12 + hdr->l_text;
            size_t refHeader = 0;
            for (int i = 0; i < hdr->n_targets; ++i) {
                char* n = hdr->target_name[i];
                refHeader += (8 + (strlen(n) + 1));
            }
            return textHeader + refHeader;
        };
        uncompressedFilePos_ = headerLength(header_.get());
    }

    void OpenGzi(size_t numThreads)
    {
        size_t actualNumThreads = numThreads;
        if (actualNumThreads == 0) {
            actualNumThreads = std::thread::hardware_concurrency();

            // if still unknown, default to single-threaded
            if (actualNumThreads == 0) actualNumThreads = 1;
        }
        gziThread_ = std::thread{&IndexedBamWriterPrivate2::RunGziThread, this, actualNumThreads};
    }

    void OpenPbi(const PbiBuilder::CompressionLevel compressionLevel, const size_t numThreads,
                 const size_t fileBufferSize)
    {
        builder_ = std::make_unique<PbiBuilder2>(bamFilename_, bamFilename_ + ".pbi",
                                                 compressionLevel, numThreads, fileBufferSize);
    }

    void RunGziThread(size_t numThreads)
    {
        //
        // This thread is the GZI index-enabled reader that trails the writer
        // thread(s). It checks for changes in the output BAM's file size &
        // reads whatever data is available. When writing is complete, it reads
        // anything that might remain & dumps the GZI index contents to disk.
        // This index is used downstream to generate records' "virtual offsets".
        //

        const auto& bamFilename = bamFilename_;
        std::unique_ptr<BGZF, HtslibBgzfDeleter> bgzf;

        struct stat st;
        int ret = 0;
        int64_t lastFileSize = 0;
        int64_t numBytesRead = 0;

        auto initBgzf = [&bgzf, &bamFilename, numThreads]() {
            bgzf.reset(bgzf_open(bamFilename.c_str(), "rb"));
            if (!bgzf) throw std::runtime_error{"could not open BAM for toy train reading"};
            bgzf_index_build_init(bgzf.get());
            if (numThreads > 1) bgzf_mt(bgzf.get(), numThreads, 256);
        };

        // main thread loop
        while (true) {
            // Quit if writer thread(s) are finished.
            if (done_) break;

            if (stat(bamFilename.c_str(), &st) != 0) {
                gziStatus_ = GziStatus::MISC_ERROR;
                return;
            }
            if (st.st_size > lastFileSize) {
                lastFileSize = st.st_size;
            } else {
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
                continue;
            }

            // Don't read unless we can guarantee we won't catch up to the end of the file.
            // Otherwise htslib will think the file has been truncated and throw errors.
            // This is a touch tricky because we're reading in multi-thread mode so htslib
            // will speculatively start grabbing blocks.  So we're going to stay *well* behind
            // This needs be made more robust.
            //
            // Note: It's worth noting that bgzf->block_clength might only be the length of the
            // compressed *payload*, meaning if there is any other header/metadata/etc on disk
            // in the actual file, our estimation of our trailing distance might be off.  If
            // this ever starts throwing exceptions we'll have to look more in to this...
            while (lastFileSize - numBytesRead > 100 * BGZF_MAX_BLOCK_SIZE) {
                // Open BAM reader if not already open.  Need to make sure we don't open it
                // until we've already established our trailing distance.
                if (!bgzf) initBgzf();

                auto result = bgzf_read_block(bgzf.get());
                if (result != 0) {
                    gziStatus_ = GziStatus::IO_ERROR;
                    return;
                }
                if (bgzf->block_length == 0) {
                    gziStatus_ = GziStatus::TRAIL_ERROR;
                    return;
                }
                numBytesRead += bgzf->block_clength;
            }

            // Only update if thigs have appreciably fallen behind
            if (lastFileSize - numBytesRead > 1.10 * maxTrailingDistance_)
                maxTrailingDistance_ = lastFileSize - numBytesRead;
        }

        // Try to open BAM if it wasn't opened in main loop.
        if (!bgzf) initBgzf();

        // Read any remaining data.
        while (true) {
            auto result = bgzf_read_block(bgzf.get());
            if (result != 0) {
                gziStatus_ = GziStatus::IO_ERROR;
                return;
            }
            if (bgzf->block_length == 0) break;
        }

        // Dump GZI contents to disk.
        const std::string gziFn{bamFilename_ + ".gzi"};
        ret = bgzf_index_dump(bgzf.get(), gziFn.c_str(), nullptr);
        if (ret != 0) gziStatus_ = GziStatus::GZI_ERROR;
    }

public:
    void Write(const BamRecord& record)
    {
// TODO: add API to auto-skip this without special compile flag
#if PBBAM_AUTOVALIDATE
        Validator::Validate(record);
#endif
        // add record & its to index builder.
        //
        // NOTE: This is the record's postiion as if it were _uncompressed_. We
        //       will return with GZI data later to transform it into BAM
        //       "virtual offset".
        //
        builder_->AddRecord(record, uncompressedFilePos_);

        const auto rawRecord = internal::BamRecordMemory::GetRawData(record);

        // update bin
        // min_shift=14 & n_lvls=5 are BAM "magic numbers"
        rawRecord->core.bin = hts_reg2bin(rawRecord->core.pos, bam_endpos(rawRecord.get()), 14, 5);

        // write record to file
        const auto ret = sam_write1(bam_.get(), header_.get(), rawRecord.get());
        if (ret <= 0) throw std::runtime_error{"could not write record"};

        // update file position
        auto recordLength = [](bam1_t* b) {
            auto* c = &b->core;

            static constexpr size_t fixedLength = 36;
            const size_t qnameLength = (c->l_qname - c->l_extranul);

            // TODO: long CIGAR handling... sigh...

            size_t remainingLength = 0;
            if (c->n_cigar <= 0xffff)
                remainingLength = (b->l_data - c->l_qname);
            else {
                const size_t cigarEnd = ((uint8_t*)bam_get_cigar(b) - b->data) + (c->n_cigar * 4);
                remainingLength = 8 + (b->l_data - cigarEnd) + 4 + (4 * c->n_cigar);
            }

            return fixedLength + qnameLength + remainingLength;
        };
        uncompressedFilePos_ += recordLength(rawRecord.get());

        // Need to handle any errors from the gzi thread, since it's not set
        // up to throw without terminating the program
        auto gstatus = gziStatus_.load();
        if (gstatus != GziStatus::GOOD) {
            if (gziStatus_.load() == GziStatus::IO_ERROR)
                throw std::runtime_error("Error in gzi thread reading from BAM file " +
                                         bamFilename_);
            if (gziStatus_.load() == GziStatus::TRAIL_ERROR)
                throw std::runtime_error(
                    "Gzi reader thread failed to properly trail when reading " + bamFilename_);
            if (gziStatus_.load() == GziStatus::GZI_ERROR)
                throw std::runtime_error("Could not dump GZI contents for indexing " +
                                         bamFilename_);
            if (gziStatus_.load() == GziStatus::MISC_ERROR)
                throw std::runtime_error("Error computing index file for " + bamFilename_);
            gziStatus_.store(GziStatus::DEAD);
        }
    }

    size_t MaxReaderLag() const { return maxTrailingDistance_; }

private:
    std::string bamFilename_;

    std::shared_ptr<bam_hdr_t> header_;
    std::unique_ptr<samFile, internal::HtslibFileDeleter> bam_;
    std::unique_ptr<PbiBuilder2> builder_;

    // used as a type of error return code for the gziThread, so
    // that errors are delayed until at least the bam file is
    // safely written to disk
    enum class GziStatus
    {
        GOOD,
        IO_ERROR,
        TRAIL_ERROR,
        GZI_ERROR,
        MISC_ERROR,
        // There was an error, but we've bubbled up the
        // information already
        DEAD
    };
    std::atomic<GziStatus> gziStatus_{GziStatus::GOOD};
    std::thread gziThread_;

    bool blockWritten_ = false;
    bool isOpen_ = false;

    std::atomic<bool> done_{false};
    std::atomic<size_t> maxTrailingDistance_{0};

    int64_t uncompressedFilePos_ = 0;
};

}  // namespace internal

IndexedBamWriter::IndexedBamWriter(const std::string& outputFilename, const BamHeader& header,
                                   const BamWriter::CompressionLevel bamCompressionLevel,
                                   const size_t numBamThreads,
                                   const PbiBuilder::CompressionLevel pbiCompressionLevel,
                                   const size_t numPbiThreads, const size_t numGziThreads,
                                   const size_t tempFileBufferSize)
    : IRecordWriter(), d_{nullptr}
{
#if PBBAM_AUTOVALIDATE
    Validator::Validate(header);
#endif
    d_ = std::make_unique<internal::IndexedBamWriterPrivate2>(
        outputFilename, internal::BamHeaderMemory::MakeRawHeader(header), bamCompressionLevel,
        numBamThreads, pbiCompressionLevel, numPbiThreads, numGziThreads, tempFileBufferSize);
}

IndexedBamWriter::~IndexedBamWriter() {}

void IndexedBamWriter::TryFlush() {}  // ignore

void IndexedBamWriter::Write(const BamRecord& record) { d_->Write(record); }

void IndexedBamWriter::Write(const BamRecordImpl& record) { d_->Write(BamRecord{record}); }

size_t IndexedBamWriter::MaxReaderLag() const { return d_->MaxReaderLag(); }

}  // namespace BAM
}  // namespace PacBio
