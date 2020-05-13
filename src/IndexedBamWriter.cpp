/// File Description
/// \file IndexedBamWriter.cpp
/// \brief Implements the IndexedBamWriter class
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/IndexedBamWriter.h"

#include <sys/stat.h>

#include <cassert>
#include <cstdint>

#include <array>
#include <atomic>
#include <condition_variable>
#include <fstream>
#include <iostream>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <tuple>
#include <type_traits>

#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <pbcopper/utility/Deleters.h>
#include <boost/numeric/conversion/cast.hpp>

#include "pbbam/BamHeader.h"
#include "pbbam/BamRecord.h"
#include "pbbam/BamRecordImpl.h"
#include "pbbam/BamWriter.h"
#include "pbbam/PbiRawData.h"
#include "pbbam/RecordType.h"
#include "pbbam/Validator.h"

#include "FileProducer.h"
#include "MemoryUtils.h"

namespace PacBio {
namespace BAM {

struct IndexedBamWriterException : public std::exception
{
    IndexedBamWriterException(std::string filename, std::string reason)
    {
        std::ostringstream s;
        s << "[pbbam] indexed BAM writer ERROR: " << reason << ":\n"
          << "  file: " << filename;
        msg_ = s.str();
    }

    const char* what() const noexcept override { return msg_.c_str(); }

    std::string msg_;
};

namespace internal {

void bgzf_write_safe2(BGZF* fp, const void* data, size_t length)
{
    const auto ret = bgzf_write(fp, data, length);
    if (ret < 0L)
        throw std::runtime_error{
            "[pbbam] indexed BAM writer ERROR: non-zero returned from bgzf_write(). Out of disk "
            "space?"};
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
            throw std::runtime_error{
                "[pbbam] indexed BAM writer ERROR: unsupported element size: " +
                std::to_string(elementSize)};
    }
}

template <typename T>
inline void WriteBgzfVector2(BGZF* fp, std::vector<T>& data)
{
    assert(fp);
    if (fp->is_be) SwapEndianness2(data);
    bgzf_write_safe2(fp, data.data(), data.size() * sizeof(T));
}

struct PbiFieldBlock2
{
    int64_t pos_;  // file position of block start
    size_t n_;     // number of entries in block
};

template <typename T>
class PbiField2
{
    constexpr static const size_t ElementSize = sizeof(T);

public:
    PbiField2(size_t maxBufferSize) : maxElementCount_{maxBufferSize / ElementSize}
    {
        buffer_.reserve(maxElementCount_);
    }

    void Add(T value) { buffer_.push_back(value); }
    bool IsFull() const { return buffer_.size() == maxElementCount_; }

    size_t maxElementCount_;
    std::vector<T> buffer_;
    std::vector<PbiFieldBlock2> blocks_;
};

class PbiReferenceDataBuilder2
{
public:
    using ReferenceRows = std::pair<int32_t, int32_t>;  // [startRow, endRow)

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
    enum class FlushMode
    {
        FORCE,
        NO_FORCE
    };

public:
    PbiBuilder2(const std::string& bamFilename, const std::string& pbiFilename,
                const PbiBuilder::CompressionLevel compressionLevel, const size_t numThreads,
                const size_t fileBufferSize)
        //                const size_t numReferenceSequences = 0
        //                const bool isCoordinateSorted = false
        : bamFilename_{bamFilename},
          pbiFilename_{pbiFilename},
          tempFilename_{pbiFilename + ".build"},
          tempFile_{std::fopen(tempFilename_.c_str(), "w+b")},
          compressionLevel_{compressionLevel},
          numThreads_{numThreads},
          rgIdField_{fileBufferSize},
          qStartField_{fileBufferSize},
          qEndField_{fileBufferSize},
          holeNumField_{fileBufferSize},
          readQualField_{fileBufferSize},
          ctxtField_{fileBufferSize},
          fileOffsetField_{fileBufferSize},
          tIdField_{fileBufferSize},
          tStartField_{fileBufferSize},
          tEndField_{fileBufferSize},
          aStartField_{fileBufferSize},
          aEndField_{fileBufferSize},
          revStrandField_{fileBufferSize},
          nMField_{fileBufferSize},
          nMMField_{fileBufferSize},
          mapQualField_{fileBufferSize},
          bcForwardField_{fileBufferSize},
          bcReverseField_{fileBufferSize},
          bcQualField_{fileBufferSize}
    {
        if (!tempFile_) throw IndexedBamWriterException{tempFilename_, "could not open temp file"};

        // TODO: setup for ref data building
    }

    void AddRecord(const BamRecord& b, const int64_t uOffset)
    {
        // ensure updated data (necessary?)
        PacBio::BAM::BamRecordMemory::UpdateRecordTags(b);
        b.ResetCachedPositions();

        // store record data & maybe flush to temp file
        AddBasicData(b, uOffset);
        AddMappedData(b);
        AddBarcodeData(b);
        AddReferenceData(b, currentRow_);
        FlushBuffers(FlushMode::NO_FORCE);

        ++currentRow_;
    }

    void AddBasicData(const BamRecord& b, const int64_t uOffset)
    {
        // read group ID
        const auto rgId = [&b]() -> int32_t {
            auto rgIdString = b.ReadGroupBaseId();
            if (rgIdString.empty()) rgIdString = MakeReadGroupId(b.MovieName(), ToString(b.Type()));
            return std::stoul(rgIdString, nullptr, 16);
        }();

        // query start/end
        const auto isCcsOrTranscript = (IsCcsOrTranscript(b.Type()));
        const int32_t qStart = (isCcsOrTranscript ? 0 : b.QueryStart());
        const int32_t qEnd = (isCcsOrTranscript ? b.Impl().SequenceLength() : b.QueryEnd());

        // add'l data
        const int32_t holeNum = (b.HasHoleNumber() ? b.HoleNumber() : 0);
        const float readAccuracy =
            (b.HasReadAccuracy() ? boost::numeric_cast<float>(b.ReadAccuracy()) : 0.0F);
        const uint8_t ctxt = (b.HasLocalContextFlags() ? b.LocalContextFlags()
                                                       : LocalContextFlags::NO_LOCAL_CONTEXT);

        // store
        rgIdField_.Add(rgId);
        qStartField_.Add(qStart);
        qEndField_.Add(qEnd);
        holeNumField_.Add(holeNum);
        ctxtField_.Add(ctxt);
        readQualField_.Add(readAccuracy);
        fileOffsetField_.Add(uOffset);
    }

    void AddMappedData(const BamRecord& b)
    {
        // alignment position
        const auto tId = b.ReferenceId();
        const auto tStart = static_cast<uint32_t>(b.ReferenceStart());
        const auto tEnd = static_cast<uint32_t>(b.ReferenceEnd());
        const auto aStart = static_cast<uint32_t>(b.AlignedStart());
        const auto aEnd = static_cast<uint32_t>(b.AlignedEnd());
        const auto isReverseStrand = [&b]() -> uint8_t {
            return (b.AlignedStrand() == Strand::REVERSE ? 1 : 0);
        }();

        // alignment quality
        const auto matchData = b.NumMatchesAndMismatches();
        const auto nM = static_cast<uint32_t>(matchData.first);
        const auto nMM = static_cast<uint32_t>(matchData.second);
        const auto mapQuality = b.MapQuality();

        if (tId >= 0) hasMappedData_ = true;

        // store
        tIdField_.Add(tId);
        tStartField_.Add(tStart);
        tEndField_.Add(tEnd);
        aStartField_.Add(aStart);
        aEndField_.Add(aEnd);
        revStrandField_.Add(isReverseStrand);
        nMField_.Add(nM);
        nMMField_.Add(nMM);
        mapQualField_.Add(mapQuality);
    }

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
        bcForwardField_.Add(bcForward);
        bcReverseField_.Add(bcReverse);
        bcQualField_.Add(bcQuality);
    }

    void AddReferenceData(const BamRecord& b, const uint32_t currentRow)
    {
        // only add if coordinate-sorted hint is set
        // update with info from refDataBuilder
        if (refDataBuilder_) {
            const auto sorted = refDataBuilder_->AddRecord(b, currentRow);
            if (!sorted) refDataBuilder_.reset();
        }
    }

    void Close()
    {
        if (isClosed_) return;

        FlushBuffers(FlushMode::FORCE);

        OpenPbiFile();
        WritePbiHeader();
        WriteFromTempFile();

        remove(tempFilename_.c_str());
        isClosed_ = true;
    }

    void OpenPbiFile()
    {
        // open file handle
        const auto mode = std::string("wb") + std::to_string(static_cast<int>(compressionLevel_));
        pbiFile_.reset(bgzf_open(pbiFilename_.c_str(), mode.c_str()));
        if (pbiFile_ == nullptr)
            throw IndexedBamWriterException{pbiFilename_, "could not open output *.pbi file"};

        // if no explicit thread count given, attempt built-in check
        size_t actualNumThreads = numThreads_;
        if (actualNumThreads == 0) {
            actualNumThreads = std::thread::hardware_concurrency();

            // if still unknown, default to single-threaded
            if (actualNumThreads == 0) actualNumThreads = 1;
        }

        // if multithreading requested, enable it
        if (actualNumThreads > 1) bgzf_mt(pbiFile_.get(), actualNumThreads, 256);
    }

    template <typename T>
    void MaybeFlushBuffer(PbiField2<T>& field, bool force)
    {
        // replace with lambda, in FlushBuffer(), once PPA can use C++14 ?
        if (field.IsFull() || force) {
            WriteToTempFile(field);
            field.buffer_.clear();
        }
    }

    void FlushBuffers(FlushMode mode)
    {
        const auto force = (mode == FlushMode::FORCE);

        MaybeFlushBuffer(rgIdField_, force);
        MaybeFlushBuffer(qStartField_, force);
        MaybeFlushBuffer(qEndField_, force);
        MaybeFlushBuffer(holeNumField_, force);
        MaybeFlushBuffer(readQualField_, force);
        MaybeFlushBuffer(ctxtField_, force);
        MaybeFlushBuffer(fileOffsetField_, force);

        MaybeFlushBuffer(tIdField_, force);
        MaybeFlushBuffer(tStartField_, force);
        MaybeFlushBuffer(tEndField_, force);
        MaybeFlushBuffer(aStartField_, force);
        MaybeFlushBuffer(aEndField_, force);
        MaybeFlushBuffer(revStrandField_, force);
        MaybeFlushBuffer(nMField_, force);
        MaybeFlushBuffer(nMMField_, force);
        MaybeFlushBuffer(mapQualField_, force);

        MaybeFlushBuffer(bcForwardField_, force);
        MaybeFlushBuffer(bcReverseField_, force);
        MaybeFlushBuffer(bcQualField_, force);
    }

    template <typename T>
    void LoadFieldBlockFromTempFile(PbiField2<T>& field, const PbiFieldBlock2& block)
    {
        // seek to block begin
        const auto ret = std::fseek(tempFile_.get(), block.pos_, SEEK_SET);
        if (ret != 0) {
            std::ostringstream s;
            s << "[pbbam] indexed BAM writer ERROR: could not seek in temp file:\n"
              << "  file: " << tempFilename_ << '\n'
              << "  offset: " << block.pos_;
            throw std::runtime_error{s.str()};
        }

        // read block elements
        field.buffer_.assign(block.n_, 0);
        const auto numElements =
            std::fread(field.buffer_.data(), sizeof(T), block.n_, tempFile_.get());

        if (numElements != block.n_) {
            std::ostringstream s;
            s << "[pbbam] indexed BAM writer ERROR: could not read expected element count:\n"
              << "  file: " << tempFilename_ << '\n'
              << "  expected: " << block.n_ << '\n'
              << "  observed: " << numElements;
            throw std::runtime_error{s.str()};
        }
    }

    template <typename T>
    void WriteField(PbiField2<T>& field)
    {
        for (const auto& block : field.blocks_) {
            LoadFieldBlockFromTempFile(field, block);
            WriteBgzfVector2(pbiFile_.get(), field.buffer_);
        }
    }

    void WriteFromTempFile()
    {
        // load from temp file, in PBI format order, and write to index

        WriteField(rgIdField_);
        WriteField(qStartField_);
        WriteField(qEndField_);
        WriteField(holeNumField_);
        WriteField(readQualField_);
        WriteField(ctxtField_);

        WriteVirtualOffsets();

        if (hasMappedData_) {
            WriteField(tIdField_);
            WriteField(tStartField_);
            WriteField(tEndField_);
            WriteField(aStartField_);
            WriteField(aEndField_);
            WriteField(revStrandField_);
            WriteField(nMField_);
            WriteField(nMMField_);
            WriteField(mapQualField_);
        }

        if (refDataBuilder_) WriteReferenceData();

        if (hasBarcodeData_) {
            WriteField(bcForwardField_);
            WriteField(bcReverseField_);
            WriteField(bcQualField_);
        }
    }

    template <typename T>
    void WriteToTempFile(PbiField2<T>& field)
    {
        if (field.buffer_.empty()) return;

        const auto pos = std::ftell(tempFile_.get());
        const auto numElements =
            std::fwrite(field.buffer_.data(), sizeof(T), field.buffer_.size(), tempFile_.get());
        field.blocks_.emplace_back(PbiFieldBlock2{pos, numElements});
    }

    void WritePbiHeader()
    {
        BGZF* bgzf = pbiFile_.get();

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

    void WriteReferenceData() { refDataBuilder_->WriteData(pbiFile_.get()); }

    std::vector<GzIndexEntry> LoadGzi()
    {
        //
        // Open GZI file & load its contents. About to use for offset transformation.
        //

        const std::string gziFn{bamFilename_ + ".gzi"};
        std::unique_ptr<FILE, Utility::FileDeleter> gziFile{fopen(gziFn.c_str(), "rb")};
        if (!gziFile) throw IndexedBamWriterException{gziFn, "could not open *.gzi file"};

        uint64_t numElements;
        if (fread(&numElements, sizeof(numElements), 1, gziFile.get()) < 1)
            throw IndexedBamWriterException{gziFn, "could not read from *.gzi file"};
        if (ed_is_big()) ed_swap_8(numElements);

        std::vector<GzIndexEntry> result;
        result.reserve(numElements);
        for (uint32_t i = 0; i < numElements; ++i) {
            GzIndexEntry entry;
            const auto vReturn = fread(&entry.vAddress, sizeof(entry.vAddress), 1, gziFile.get());
            const auto uReturn = fread(&entry.uAddress, sizeof(entry.uAddress), 1, gziFile.get());
            if (vReturn < 1 || uReturn < 1)
                throw IndexedBamWriterException{gziFn, "could not read from *.gzi file"};

            if (ed_is_big()) {
                ed_swap_8(entry.vAddress);
                ed_swap_8(entry.uAddress);
            }
            result.push_back(std::move(entry));
        }
        return result;
    }

    void WriteVirtualOffsets()
    {
        auto index = LoadGzi();
        if (index.empty())
            throw IndexedBamWriterException{bamFilename_ + ".gzi", "empty *.gzi file contents"};
        std::sort(index.begin(), index.end(),
                  [](const GzIndexEntry& lhs, const GzIndexEntry& rhs) -> bool {
                      return lhs.uAddress < rhs.uAddress;
                  });

        size_t k = 0;
        for (const auto& block : fileOffsetField_.blocks_) {
            LoadFieldBlockFromTempFile(fileOffsetField_, block);

            // transform offsets from GZI
            for (size_t j = 0; j < fileOffsetField_.buffer_.size(); ++j) {
                while ((k < index.size() - 1) && (static_cast<uint64_t>(index.at(k + 1).uAddress) <=
                                                  fileOffsetField_.buffer_[j])) {
                    ++k;
                }
                const GzIndexEntry& e = index.at(k);
                const int64_t uOffset = fileOffsetField_.buffer_[j] - e.uAddress;
                const auto result = ((e.vAddress << 16) | uOffset);
                fileOffsetField_.buffer_[j] = result;
            }
            WriteBgzfVector2(pbiFile_.get(), fileOffsetField_.buffer_);
        }
    }

private:
    // file info
    std::string bamFilename_;
    std::string pbiFilename_;
    std::string tempFilename_;
    std::unique_ptr<FILE, Utility::FileDeleter> tempFile_;
    std::unique_ptr<BGZF, HtslibBgzfDeleter> pbiFile_;
    PbiBuilder::CompressionLevel compressionLevel_;
    size_t numThreads_;

    // PBI field buffers
    PbiField2<int32_t> rgIdField_;
    PbiField2<int32_t> qStartField_;
    PbiField2<int32_t> qEndField_;
    PbiField2<int32_t> holeNumField_;
    PbiField2<float> readQualField_;
    PbiField2<uint8_t> ctxtField_;
    PbiField2<uint64_t> fileOffsetField_;
    PbiField2<int32_t> tIdField_;
    PbiField2<uint32_t> tStartField_;
    PbiField2<uint32_t> tEndField_;
    PbiField2<uint32_t> aStartField_;
    PbiField2<uint32_t> aEndField_;
    PbiField2<uint8_t> revStrandField_;
    PbiField2<uint32_t> nMField_;
    PbiField2<uint32_t> nMMField_;
    PbiField2<uint8_t> mapQualField_;
    PbiField2<int16_t> bcForwardField_;
    PbiField2<int16_t> bcReverseField_;
    PbiField2<int8_t> bcQualField_;

    // reference data
    std::unique_ptr<PbiReferenceDataBuilder2> refDataBuilder_;

    // tracking data
    uint32_t currentRow_ = 0;
    bool isClosed_ = false;
    bool hasBarcodeData_ = false;
    bool hasMappedData_ = false;
};

}  // namespace internal

class IndexedBamWriter::IndexedBamWriterPrivate2  //: public internal::FileProducer
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

        remove(std::string{bamFilename_ + ".gzi"}.c_str());
        isOpen_ = false;
    }

    void CloseBam()
    {
        const auto ret = bgzf_flush(bam_.get()->fp.bgzf);
        std::ignore = ret;
        bam_.reset();
    }

    void CloseGzi()
    {
        done_ = true;
        gziThread_.join();

        // TODO: remove GZI file, leaving now for debubging
    }

    void ClosePbi() { builder_->Close(); }

    void OpenBam(const BamWriter::CompressionLevel compressionLevel, const size_t numThreads)
    {
        //
        // TODO: Compression level & numThreads are hardcoded here. Ok for
        //       prototyping but need to be tune-able via API.
        //

        if (!header_) throw IndexedBamWriterException{bamFilename_, "null header provided"};

        // open output BAM
        const auto usingFilename = bamFilename_;
        const auto mode = std::string("wb") + std::to_string(static_cast<int>(compressionLevel));
        bam_.reset(sam_open(usingFilename.c_str(), mode.c_str()));
        if (!bam_)
            throw IndexedBamWriterException{usingFilename, "could not open file for writing"};

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
        if (ret != 0) throw IndexedBamWriterException{usingFilename, "could not write header"};
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
        builder_ = std::make_unique<internal::PbiBuilder2>(
            bamFilename_, bamFilename_ + ".pbi", compressionLevel, numThreads, fileBufferSize);
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
            if (!bgzf)
                throw IndexedBamWriterException{bamFilename, "could not open trailing BAM reader"};
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

        const auto rawRecord = BamRecordMemory::GetRawData(record);

        // update bin
        // min_shift=14 & n_lvls=5 are BAM "magic numbers"
        rawRecord->core.bin = hts_reg2bin(rawRecord->core.pos, bam_endpos(rawRecord.get()), 14, 5);

        // write record to file
        const auto ret = sam_write1(bam_.get(), header_.get(), rawRecord.get());
        if (ret <= 0) throw IndexedBamWriterException{bamFilename_, "could not write record"};

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
                throw IndexedBamWriterException{bamFilename_,
                                                "error in gzi thread reading from BAM file"};
            if (gziStatus_.load() == GziStatus::TRAIL_ERROR)
                throw IndexedBamWriterException{
                    bamFilename_, "gzi reader thread failed to properly trail when reading"};
            if (gziStatus_.load() == GziStatus::GZI_ERROR)
                throw IndexedBamWriterException{bamFilename_,
                                                "could not dump GZI contents for indexing"};
            if (gziStatus_.load() == GziStatus::MISC_ERROR)
                throw IndexedBamWriterException{bamFilename_, "error computing index file"};
            gziStatus_.store(GziStatus::DEAD);
        }
    }

    size_t MaxReaderLag() const { return maxTrailingDistance_; }

private:
    std::string bamFilename_;

    std::shared_ptr<bam_hdr_t> header_;
    std::unique_ptr<samFile, HtslibFileDeleter> bam_;
    std::unique_ptr<internal::PbiBuilder2> builder_;

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

    bool isOpen_ = false;

    std::atomic<bool> done_{false};
    std::atomic<size_t> maxTrailingDistance_{0};

    int64_t uncompressedFilePos_ = 0;
};

static_assert(!std::is_copy_constructible<IndexedBamWriter>::value,
              "IndexedBamWriter(const IndexedBamWriter&) is not = delete");
static_assert(!std::is_copy_assignable<IndexedBamWriter>::value,
              "IndexedBamWriter& operator=(const IndexedBamWriter&) is not = delete");

IndexedBamWriter::IndexedBamWriter(const std::string& outputFilename, const BamHeader& header,
                                   const BamWriter::CompressionLevel bamCompressionLevel,
                                   const size_t numBamThreads,
                                   const PbiBuilder::CompressionLevel pbiCompressionLevel,
                                   const size_t numPbiThreads, const size_t numGziThreads,
                                   const size_t tempFileBufferSize)
    : IRecordWriter(), d_{nullptr}
{
    if (tempFileBufferSize % 8 != 0)
        throw std::runtime_error{
            "[pbbam] indexed BAM writer ERROR: invalid buffer size for PBI builder (" +
            std::to_string(tempFileBufferSize) + "). Must be a multiple of 8."};

#if PBBAM_AUTOVALIDATE
    Validator::Validate(header);
#endif
    d_ = std::make_unique<IndexedBamWriterPrivate2>(
        outputFilename, BamHeaderMemory::MakeRawHeader(header), bamCompressionLevel, numBamThreads,
        pbiCompressionLevel, numPbiThreads, numGziThreads, tempFileBufferSize);
}

IndexedBamWriter::IndexedBamWriter(IndexedBamWriter&&) noexcept = default;

IndexedBamWriter& IndexedBamWriter::operator=(IndexedBamWriter&&) noexcept = default;

IndexedBamWriter::~IndexedBamWriter() = default;

void IndexedBamWriter::Write(const BamRecord& record) { d_->Write(record); }

void IndexedBamWriter::Write(const BamRecordImpl& record) { d_->Write(BamRecord{record}); }

size_t IndexedBamWriter::MaxReaderLag() const { return d_->MaxReaderLag(); }

}  // namespace BAM
}  // namespace PacBio
