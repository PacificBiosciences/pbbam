// File Description
/// \file CCSPbiBuilder.cpp
/// \brief Implements the CCSPbiBuilder.cpp class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/ccs/CCSPbiBuilder.h"

#include <cstdio>
#include <cstdlib>

#include <stdexcept>
#include <thread>
#include <vector>

#include <htslib/bgzf.h>

#include "pbbam/PbiBuilder.h"
#include "pbbam/PbiFile.h"
#include "pbbam/ccs/CCSHeader.h"
#include "pbbam/ccs/CCSRecord.h"

#include "MemoryUtils.h"

namespace PacBio {
namespace CCS {
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
            throw std::runtime_error{"CCSPbiBuilder: unsupported element size (" +
                                     std::to_string(elementSize) + ")"};
    }
}

void bgzf_write_safe(BGZF* fp, const void* data, size_t length)
{
    const auto ret = bgzf_write(fp, data, length);
    if (ret < 0L)
        throw std::runtime_error{
            "CCSPbiBuilder: non-zero returned from bgzf_write(). Out of disk space?"};
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

}  // namespace internal

class CCSPbiBuilder::CCSPbiBuilderPrivate
{
    enum class FlushMode
    {
        FORCE,
        NO_FORCE
    };

    // TODO: Make this tweak-able, a la IndexedBamWriter's buffers
    constexpr static const size_t MaxBufferSize = 0x10000;

public:
    CCSPbiBuilderPrivate(const std::string& pbiFilename, const std::string& movieName,
                         const CCSPbiBuilderConfig& config)
        : pbiFilename_{pbiFilename}
        , tempFilename_{pbiFilename + ".build"}
        , tempFile_{std::fopen(tempFilename_.c_str(), "w+b")}
        , compressionLevel_{config.CompressionLevel}
        , numThreads_{config.NumThreads}
        , rgIdField_{MaxBufferSize}
        , qStartField_{MaxBufferSize}
        , qEndField_{MaxBufferSize}
        , holeNumField_{MaxBufferSize}
        , readQualField_{MaxBufferSize}
        , ctxtField_{MaxBufferSize}
        , fileOffsetField_{MaxBufferSize}
    {
        movieName_ = movieName;
        rgId_ = BAM::ReadGroupInfo::IdToInt(BAM::MakeReadGroupId(movieName, "SUBREAD"));
    }

    void AddRecord(const CCSRecord& record)
    {
        rgIdField_.Add(rgId_);
        qStartField_.Add(record.QueryStart);
        qEndField_.Add(record.QueryEnd);
        holeNumField_.Add(record.HoleNumber);
        ctxtField_.Add(record.LocalContextFlags);
        readQualField_.Add(record.Accuracy);
        fileOffsetField_.Add(-1);

        FlushBuffers(FlushMode::NO_FORCE);
        ++currentRow_;
    }

    void Close()
    {
        if (isClosed_) return;

        FlushBuffers(FlushMode::FORCE);

        OpenPbiFile();
        WritePbiHeader();
        WriteFromTempFile();

        std::remove(tempFilename_.c_str());
        isClosed_ = true;
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

    void OpenPbiFile()
    {
        // open file handle
        const auto mode = std::string("wb") + std::to_string(static_cast<int>(compressionLevel_));
        pbiFile_.reset(bgzf_open(pbiFilename_.c_str(), mode.c_str()));
        if (pbiFile_ == nullptr)
            throw std::runtime_error{"CCSPbiBuilder: could not open file for writing: " +
                                     pbiFilename_};

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

    void WritePbiHeader()
    {
        BGZF* bgzf = pbiFile_.get();

        // 'magic' string
        static constexpr const std::array<char, 4> magic{{'P', 'B', 'I', '\1'}};
        internal::bgzf_write_safe(bgzf, magic.data(), 4);

        PacBio::BAM::PbiFile::Sections sections = PacBio::BAM::PbiFile::BASIC;
        // version, pbi_flags, & n_reads
        auto version = static_cast<uint32_t>(PacBio::BAM::PbiFile::CurrentVersion);
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

    template <typename T>
    void LoadFieldBlockFromTempFile(internal::PbiField<T>& field,
                                    const internal::PbiFieldBlock& block)
    {
        // seek to block begin
        const auto ret = std::fseek(tempFile_.get(), block.pos_, SEEK_SET);
        if (ret != 0)
            throw std::runtime_error{"CCSPbiBuilder: could not seek in temp file: " +
                                     tempFilename_ + ", offset: " + std::to_string(block.pos_)};

        // read block elements
        field.buffer_.assign(block.n_, 0);
        const auto numElements =
            std::fread(field.buffer_.data(), sizeof(T), block.n_, tempFile_.get());

        if (numElements != block.n_)
            throw std::runtime_error{
                "CCSPbiBuilder: could not read element count from temp file: " + tempFilename_};
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
    }

    // file info
    std::string bamFilename_;
    std::string pbiFilename_;
    std::string tempFilename_;
    std::unique_ptr<FILE, PacBio::BAM::FileDeleter> tempFile_;
    std::unique_ptr<BGZF, PacBio::BAM::HtslibBgzfDeleter> pbiFile_;
    PacBio::BAM::PbiBuilder::CompressionLevel compressionLevel_ =
        PacBio::BAM::PbiBuilder::DefaultCompression;
    size_t numThreads_;

    // PBI field buffers
    internal::PbiField<int32_t> rgIdField_;
    internal::PbiField<int32_t> qStartField_;
    internal::PbiField<int32_t> qEndField_;
    internal::PbiField<int32_t> holeNumField_;
    internal::PbiField<float> readQualField_;
    internal::PbiField<uint8_t> ctxtField_;
    internal::PbiField<uint64_t> fileOffsetField_;

    std::string movieName_;
    int32_t rgId_;
    uint32_t currentRow_ = 0;
    bool isClosed_ = false;
};

CCSPbiBuilder::CCSPbiBuilder(const std::string& pbiFilename, const std::string& movieName,
                             const CCSPbiBuilderConfig& config)
    : d_{std::make_unique<CCSPbiBuilderPrivate>(pbiFilename, movieName, config)}
{
}

CCSPbiBuilder::CCSPbiBuilder(const std::string& pbiFilename, const CCSHeader& header,
                             const CCSPbiBuilderConfig& config)
    : CCSPbiBuilder{pbiFilename, header.MovieName, config}
{
}

CCSPbiBuilder::~CCSPbiBuilder() = default;

void CCSPbiBuilder::AddRecord(const CCSRecord& record) { d_->AddRecord(record); }

void CCSPbiBuilder::Close() { d_->Close(); }

const std::string& CCSPbiBuilder::MovieName() const { return d_->movieName_; }

}  // namespace CCS
}  // namespace PacBio