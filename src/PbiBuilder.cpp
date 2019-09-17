// File Description
/// \file PbiBuilder.cpp
/// \brief Implements the PbiBuilder class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/PbiBuilder.h"

#include <cstddef>
#include <cstdint>
#include <cstdio>

#include <array>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <tuple>

#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>

#include <boost/numeric/conversion/cast.hpp>

#include <pbcopper/utility/Deleters.h>

#include "pbbam/BamRecord.h"
#include "pbbam/BamRecordImpl.h"
#include "pbbam/PbiRawData.h"
#include "pbbam/RecordType.h"

#include "MemoryUtils.h"

namespace PacBio {
namespace BAM {

struct PbiBuilderException : public std::exception
{
    PbiBuilderException(std::string filename, std::string reason)
    {
        std::ostringstream s;
        s << "[pbbam] PBI index builder ERROR: " << reason << ":\n"
          << "  file: " << filename;
        msg_ = s.str();
    }

    const char* what() const noexcept override { return msg_.c_str(); }

    std::string msg_;
};

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
            throw std::runtime_error{"[pbbam] PBI index builder ERROR: unsupported element size (" +
                                     std::to_string(elementSize) + ")"};
    }
}

void bgzf_write_safe(BGZF* fp, const void* data, size_t length)
{
    const auto ret = bgzf_write(fp, data, length);
    if (ret < 0L)
        throw std::runtime_error{
            "[pbbam] PBI index builder ERROR: non-zero returned from bgzf_write(). Out of disk "
            "space?"};
}

template <typename T>
inline void WriteBgzfVector(BGZF* fp, std::vector<T>& data)
{
    assert(fp);
    if (fp->is_be) SwapEndianness(data);
    bgzf_write_safe(fp, &data[0], data.size() * sizeof(T));
}

struct PbiFieldBlock
{
    int64_t pos_;  // file position of block start
    size_t n_;     // number of entries in block
};

template <typename T>
class PbiField
{
    constexpr static const size_t ElementSize = sizeof(T);

public:
    PbiField(size_t maxBufferSize) : maxElementCount_{maxBufferSize / ElementSize}
    {
        buffer_.reserve(maxElementCount_);
    }

    void Add(T value) { buffer_.push_back(value); }
    bool IsFull() const { return buffer_.size() == maxElementCount_; }

    size_t maxElementCount_;
    std::vector<T> buffer_;
    std::vector<PbiFieldBlock> blocks_;
};
// --------------------------
// PbiReferenceDataBuilder
// --------------------------

class PbiReferenceDataBuilder
{
public:
    using ReferenceRows = std::pair<int32_t, int32_t>;  // [startRow, endRow)

    explicit PbiReferenceDataBuilder(const size_t numReferenceSequences);

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
    rawReferenceEntries_[PbiReferenceEntry::UNMAPPED_ID] = PbiReferenceEntry{};
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
            PbiReferenceEntry& currentEntry = rawReferenceEntries_.at(static_cast<uint32_t>(tId));
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

PbiRawReferenceData PbiReferenceDataBuilder::Result() const
{
    // PbiReferenceEntries will be sorted thanks to std::map
    // tId will be at end since we're sorting on the uint cast of -1
    PbiRawReferenceData result;
    result.entries_.reserve(rawReferenceEntries_.size());
    for (const auto& entry : rawReferenceEntries_)
        result.entries_.push_back(entry.second);
    return result;
}

void PbiReferenceDataBuilder::WriteData(BGZF* bgzf)
{
    const auto refData = Result();

    // num_refs
    uint32_t numRefs = refData.entries_.size();
    if (bgzf->is_be) numRefs = ed_swap_4(numRefs);
    internal::bgzf_write_safe(bgzf, &numRefs, 4);

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
        internal::bgzf_write_safe(bgzf, &tId, 4);
        internal::bgzf_write_safe(bgzf, &beginRow, 4);
        internal::bgzf_write_safe(bgzf, &endRow, 4);
    }
}

}  // namespace internal

// --------------------------------------------
// PbiBuilderPrivate - builder implementation
// --------------------------------------------

// TODO: Come back to refseqs, sorting, etc

// TODO: We **NEED** to sync this up with the builder in IndexedBamWriter. They
//       differ slightly but should be shareable.

class PbiBuilder::PbiBuilderPrivate
{
    enum class FlushMode
    {
        FORCE,
        NO_FORCE
    };

    // TODO: Make this tweak-able, a la IndexedBamWriter's buffers
    constexpr static const size_t MaxBufferSize = 0x10000;

public:
    PbiBuilderPrivate(const std::string& pbiFilename, const size_t numReferenceSequences,
                      const bool isCoordinateSorted,
                      const PbiBuilder::CompressionLevel compressionLevel, const size_t numThreads)
        : pbiFilename_{pbiFilename}
        , tempFilename_{pbiFilename + ".build"}
        , tempFile_{std::fopen(tempFilename_.c_str(), "w+b")}
        , compressionLevel_{compressionLevel}
        , numThreads_{numThreads}
        , rgIdField_{MaxBufferSize}
        , qStartField_{MaxBufferSize}
        , qEndField_{MaxBufferSize}
        , holeNumField_{MaxBufferSize}
        , readQualField_{MaxBufferSize}
        , ctxtField_{MaxBufferSize}
        , fileOffsetField_{MaxBufferSize}
        , tIdField_{MaxBufferSize}
        , tStartField_{MaxBufferSize}
        , tEndField_{MaxBufferSize}
        , aStartField_{MaxBufferSize}
        , aEndField_{MaxBufferSize}
        , revStrandField_{MaxBufferSize}
        , nMField_{MaxBufferSize}
        , nMMField_{MaxBufferSize}
        , mapQualField_{MaxBufferSize}
        , bcForwardField_{MaxBufferSize}
        , bcReverseField_{MaxBufferSize}
        , bcQualField_{MaxBufferSize}
    {
        if (!tempFile_) throw PbiBuilderException{tempFilename_, "could not open temp file"};
        if (isCoordinateSorted && numReferenceSequences > 0)
            refDataBuilder_ =
                std::make_unique<internal::PbiReferenceDataBuilder>(numReferenceSequences);
    }

    ~PbiBuilderPrivate() noexcept
    {
        if (!isClosed_) {
            try {
                Close();
            } catch (...) {
                // swallow any exceptions & remain no-throw from dtor
            }
        }
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
            throw PbiBuilderException{pbiFilename_, "could not open file for writing"};

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
    void MaybeFlushBuffer(internal::PbiField<T>& field, bool force)
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
    void LoadFieldBlockFromTempFile(internal::PbiField<T>& field,
                                    const internal::PbiFieldBlock& block)
    {
        // seek to block begin
        const auto ret = std::fseek(tempFile_.get(), block.pos_, SEEK_SET);
        if (ret != 0) {
            std::ostringstream s;
            s << "[pbbam] PBI index builder ERROR: could not seek in temp file:\n"
              << "  file: " << tempFilename_ << '\n'
              << "  offset: " << block.pos_;
            throw std::runtime_error{s.str()};
        }

        // read block elements
        field.buffer_.assign(block.n_, 0);
        const auto numElements =
            std::fread(field.buffer_.data(), sizeof(T), block.n_, tempFile_.get());

        if (numElements != block.n_)
            throw PbiBuilderException{tempFilename_, "could not read element count from temp file"};
    }

    template <typename T>
    void WriteField(internal::PbiField<T>& field)
    {
        for (const auto& block : field.blocks_) {
            LoadFieldBlockFromTempFile(field, block);
            internal::WriteBgzfVector(pbiFile_.get(), field.buffer_);
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
        WriteField(fileOffsetField_);

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
    void WriteToTempFile(internal::PbiField<T>& field)
    {
        if (field.buffer_.empty()) return;

        const auto pos = std::ftell(tempFile_.get());
        const auto numElements =
            std::fwrite(field.buffer_.data(), sizeof(T), field.buffer_.size(), tempFile_.get());
        field.blocks_.emplace_back(internal::PbiFieldBlock{pos, numElements});
    }

    void WritePbiHeader()
    {
        BGZF* bgzf = pbiFile_.get();

        // 'magic' string
        static constexpr const std::array<char, 4> magic{{'P', 'B', 'I', '\1'}};
        internal::bgzf_write_safe(bgzf, magic.data(), 4);

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
        internal::bgzf_write_safe(bgzf, &version, 4);
        internal::bgzf_write_safe(bgzf, &pbi_flags, 2);
        internal::bgzf_write_safe(bgzf, &numReads, 4);

        // reserved space
        char reserved[18];
        memset(reserved, 0, 18);
        internal::bgzf_write_safe(bgzf, reserved, 18);
    }

    void WriteReferenceData() { refDataBuilder_->WriteData(pbiFile_.get()); }

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
    internal::PbiField<int32_t> rgIdField_;
    internal::PbiField<int32_t> qStartField_;
    internal::PbiField<int32_t> qEndField_;
    internal::PbiField<int32_t> holeNumField_;
    internal::PbiField<float> readQualField_;
    internal::PbiField<uint8_t> ctxtField_;
    internal::PbiField<uint64_t> fileOffsetField_;
    internal::PbiField<int32_t> tIdField_;
    internal::PbiField<uint32_t> tStartField_;
    internal::PbiField<uint32_t> tEndField_;
    internal::PbiField<uint32_t> aStartField_;
    internal::PbiField<uint32_t> aEndField_;
    internal::PbiField<uint8_t> revStrandField_;
    internal::PbiField<uint32_t> nMField_;
    internal::PbiField<uint32_t> nMMField_;
    internal::PbiField<uint8_t> mapQualField_;
    internal::PbiField<int16_t> bcForwardField_;
    internal::PbiField<int16_t> bcReverseField_;
    internal::PbiField<int8_t> bcQualField_;

    // reference data
    std::unique_ptr<internal::PbiReferenceDataBuilder> refDataBuilder_;

    // tracking data
    uint32_t currentRow_ = 0;
    bool isClosed_ = false;
    bool hasBarcodeData_ = false;
    bool hasMappedData_ = false;
};

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
    : d_{std::make_unique<PbiBuilderPrivate>(pbiFilename, numReferenceSequences, isCoordinateSorted,
                                             compressionLevel, numThreads)}
{
}

PbiBuilder::~PbiBuilder() noexcept = default;

void PbiBuilder::AddRecord(const BamRecord& record, const int64_t vOffset)
{
    d_->AddRecord(record, vOffset);
}

void PbiBuilder::Close() { d_->Close(); }

}  // namespace BAM
}  // namespace PacBio
