#ifndef PBBAM_PBIBUILDERBBASE_H
#define PBBAM_PBIBUILDERBBASE_H

#include "PbbamInternalConfig.h"

#include <pbbam/BamRecord.h>
#include <pbbam/PbiRawData.h>
#include "ErrnoReason.h"
#include "FileProducer.h"
#include "MemoryUtils.h"

#include <pbcopper/data/Position.h>
#include <pbcopper/utility/Deleters.h>

#include <boost/numeric/conversion/cast.hpp>

#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>

#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <type_traits>
#include <vector>

#include <cassert>
#include <cctype>
#include <cstddef>

namespace PacBio {
namespace BAM {

enum class FlushMode
{
    FORCE,
    NO_FORCE
};

struct PbiBuilderException : public std::exception
{
    PbiBuilderException(std::string filename, std::string reason)
    {
        std::ostringstream s;
        s << "[pbbam] PBI index builder ERROR: " << reason << ":\n"
          << "  file: " << filename;
        MaybePrintErrnoReason(s);
        msg_ = s.str();
    }

    const char* what() const noexcept override { return msg_.c_str(); }

    std::string msg_;
};

struct IndexedBamWriterException : public std::exception
{
    IndexedBamWriterException(std::string filename, std::string reason)
    {
        std::ostringstream s;
        s << "[pbbam] indexed BAM writer ERROR: " << reason << ":\n"
          << "  file: " << filename;
        MaybePrintErrnoReason(s);
        msg_ = s.str();
    }

    const char* what() const noexcept override { return msg_.c_str(); }

    std::string msg_;
};

template <typename T>
void SwapEndianness(std::vector<T>& data)
{
    const std::size_t elementSize = sizeof(T);
    const std::size_t numReads = data.size();
    switch (elementSize) {
        case 1:
            break;  // no swapping necessary
        case 2:
            for (std::size_t i = 0; i < numReads; ++i) {
                ed_swap_2p(&data[i]);
            }
            break;
        case 4:
            for (std::size_t i = 0; i < numReads; ++i) {
                ed_swap_4p(&data[i]);
            }
            break;
        case 8:
            for (std::size_t i = 0; i < numReads; ++i) {
                ed_swap_8p(&data[i]);
            }
            break;
        default:
            throw std::runtime_error{"[pbbam] PBI index builder ERROR: unsupported element size (" +
                                     std::to_string(elementSize) + ")"};
    }
}

inline void bgzf_write_safe(BGZF* fp, const void* data, std::size_t length)
{
    const auto ret = bgzf_write(fp, data, length);
    if (ret < 0L) {
        std::ostringstream msg;
        msg << "[pbbam] PBI index builder ERROR: could not write to BGZF file";
        MaybePrintErrnoReason(msg);
        throw std::runtime_error{msg.str()};
    }
}

template <typename T>
void WriteBgzfVector(BGZF* fp, std::vector<T>& data)
{
    assert(fp);
    if (fp->is_be) {
        SwapEndianness(data);
    }
    bgzf_write_safe(fp, &data[0], data.size() * sizeof(T));
}

struct PbiFieldBlock
{
    int64_t pos_;    // file position of block start
    std::size_t n_;  // number of entries in block
};

template <typename T>
class PbiField
{
    static constexpr std::size_t ELEMENT_SIZE = sizeof(T);
    static_assert(ELEMENT_SIZE > 0);

public:
    PbiField(std::size_t maxBufferSize = 0) : maxElementCount_{maxBufferSize / ELEMENT_SIZE}
    {
        buffer_.reserve(maxElementCount_);
    }

    void Add(T value) { buffer_.push_back(value); }
    bool IsFull() const { return buffer_.size() == maxElementCount_; }

    std::size_t maxElementCount_;
    std::vector<T> buffer_;
    std::vector<PbiFieldBlock> blocks_;
};

class PbiReferenceDataBuilder
{
public:
    using ReferenceRows = std::pair<int32_t, int32_t>;  // [startRow, endRow)

    explicit PbiReferenceDataBuilder(std::size_t numReferenceSequences);

    bool AddRecord(const BamRecord& record, int32_t rowNumber);
    PbiRawReferenceData Result() const;
    void WriteData(BGZF* bgzf);

private:
    int32_t lastRefId_ = -1;
    Data::Position lastPos_ = -1;
    std::map<uint32_t, PbiReferenceEntry> rawReferenceEntries_;
};

struct PbiBuilderBase
{
    PbiBuilderBase() = delete;
    explicit PbiBuilderBase(const std::string& pbiFilename,
                            PbiBuilder::CompressionLevel compressionLevel, std::size_t numThreads,
                            std::size_t bufferSize);
    virtual ~PbiBuilderBase() noexcept;

    void AddBarcodeData(const BamRecord& b);
    void AddBasicData(const BamRecord& b, int64_t uOffset);
    void AddMappedData(const BamRecord& b);
    void AddRecord(const BamRecord& b, int64_t uOffset);
    void AddReferenceData(const BamRecord& b, uint32_t currentRow);
    void Close();
    void FlushBuffers(FlushMode mode);
    void OpenPbiFile();
    void WriteFromTempFile();
    void WritePbiHeader();
    void WriteReferenceData();

    template <typename T>
    void LoadFieldBlockFromTempFile(PbiField<T>& field, const PbiFieldBlock& block)
    {
        // seek to block begin
        const auto ret = std::fseek(tempFile_.get(), block.pos_, SEEK_SET);
        if (ret != 0) {
            std::ostringstream s;
            s << "[pbbam] PBI index builder ERROR: could not seek in temp file:\n"
              << "  file: " << tempFilename_ << '\n'
              << "  offset: " << block.pos_;
            MaybePrintErrnoReason(s);
            throw std::runtime_error{s.str()};
        }

        // read block elements
        field.buffer_.assign(block.n_, 0);
        const auto numElements =
            std::fread(field.buffer_.data(), sizeof(T), block.n_, tempFile_.get());

        if (numElements != block.n_) {
            std::ostringstream msg;
            msg << "[pbbam] PBI builder ERROR: could not read element count from temp file\n"
                << "  expected: " << block.n_ << '\n'
                << "  observed: " << numElements << '\n'
                << "  file: " << tempFilename_ << '\n';
            MaybePrintErrnoReason(msg);
            throw std::runtime_error{msg.str()};
        }
    }

    template <typename T>
    void MaybeFlushBuffer(PbiField<T>& field, bool force)
    {
        // replace with lambda, in FlushBuffer(), once PPA can use C++14 ?
        if (field.IsFull() || force) {
            WriteToTempFile(field);
            field.buffer_.clear();
        }
    }

    template <typename T>
    void WriteField(PbiField<T>& field)
    {
        for (const auto& block : field.blocks_) {
            LoadFieldBlockFromTempFile(field, block);
            WriteBgzfVector(pbiFile_.get(), field.buffer_);
        }
    }

    template <typename T>
    void WriteToTempFile(PbiField<T>& field)
    {
        if (field.buffer_.empty()) {
            return;
        }

        const auto pos = std::ftell(tempFile_.get());
        const auto numElements =
            std::fwrite(field.buffer_.data(), sizeof(T), field.buffer_.size(), tempFile_.get());
        field.blocks_.emplace_back(PbiFieldBlock{pos, numElements});
    }

    virtual void WriteVirtualOffsets() { WriteField(fileOffsetField_); }

    // file/general info
    std::string pbiFilename_;
    std::string tempFilename_;
    std::unique_ptr<FILE, Utility::FileDeleter> tempFile_;
    std::unique_ptr<BGZF, HtslibBgzfDeleter> pbiFile_;
    PbiBuilder::CompressionLevel compressionLevel_;
    std::size_t numThreads_;

    // PBI field buffers
    PbiField<int32_t> rgIdField_;
    PbiField<int32_t> qStartField_;
    PbiField<int32_t> qEndField_;
    PbiField<int32_t> holeNumField_;
    PbiField<float> readQualField_;
    PbiField<uint8_t> ctxtField_;
    PbiField<uint64_t> fileOffsetField_;
    PbiField<int32_t> tIdField_;
    PbiField<uint32_t> tStartField_;
    PbiField<uint32_t> tEndField_;
    PbiField<uint32_t> aStartField_;
    PbiField<uint32_t> aEndField_;
    PbiField<uint8_t> revStrandField_;
    PbiField<uint32_t> nMField_;
    PbiField<uint32_t> nMMField_;
    PbiField<uint8_t> mapQualField_;
    PbiField<uint32_t> nInsOpsField_;
    PbiField<uint32_t> nDelOpsField_;
    PbiField<int16_t> bcForwardField_;
    PbiField<int16_t> bcReverseField_;
    PbiField<int8_t> bcQualField_;

    // reference data
    std::unique_ptr<PbiReferenceDataBuilder> refDataBuilder_;

    // tracking data
    uint32_t currentRow_ = 0;
    bool isClosed_ = false;
    bool hasBarcodeData_ = false;
    bool hasMappedData_ = false;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_PBIBUILDERBBASE_H
