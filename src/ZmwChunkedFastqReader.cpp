#include "PbbamInternalConfig.h"

#include <pbbam/ZmwChunkedFastqReader.h>

#include <cassert>
#include <cstdio>

#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include <htslib/kseq.h>

#include <pbbam/FaiIndex.h>
#include <pbbam/FormatUtils.h>

#include "MemoryUtils.h"
#include "ZmwChunkedFastxBgzfReader.h"
#include "ZmwChunkedFastxReaderImpl.h"
#include "ZmwChunkedFastxTextReader.h"

namespace PacBio {
namespace BAM {
namespace {

std::unique_ptr<ZmwChunkedFastxReaderImpl> MakeFastqReaderImpl(std::string filename,
                                                               const size_t numChunks)
{
    // validate extension
    if (!FormatUtils::IsFastqFilename(filename)) {
        std::ostringstream msg;
        msg << "[pbbam] chunked FASTQ reader ERROR: not a recognized FASTQ extension:\n"
            << "  file: " << filename;
        throw std::runtime_error{msg.str()};
    }

    // determine subsequence "loader" from compression type: plain-text, bgzf, or unsupported
    const auto compressionType = FormatUtils::CompressionType(filename);
    switch (compressionType) {

        case HtslibCompression::NONE:
            return std::make_unique<ZmwChunkedFastxTextReader>(std::move(filename), numChunks);
        case HtslibCompression::BGZIP:
            return std::make_unique<ZmwChunkedFastxBgzfReader>(std::move(filename), numChunks);

        case HtslibCompression::GZIP: {
            std::ostringstream msg;
            msg << "[pbbam] chunked FASTQ reader ERROR: random-access is not supported in standard "
                   "gzip format\n"
                << "  file: " << filename << "\n\n"
                << "Compressed files must be bgzipped, with accompanying *.gzi index.\n\n"
                << "To keep the original gzipped file unchanged:\n"
                << "  $ gunzip -c " << filename << " > <unzipped_file>\n"
                << "or discard the gzipped file:\n"
                << "  $ gunzip " << filename << '\n'
                << '\n'
                << "Re-compress & create *.gzi index:\n"
                << "  $ bgzip --index <unzipped_file>\n\n";
            throw std::runtime_error{msg.str()};
        }
        default:
            assert(false);  // should never get here, the way htslib currently determines type
            std::ostringstream msg;
            msg << "[pbbam] chunked FASTQ reader ERROR: could not determine compression type:\n"
                << "  file: " << filename;
            throw std::runtime_error{msg.str()};
    }
}

}  // namespace

class ZmwChunkedFastqReader::ZmwChunkedFastqReaderPrivate
{
public:
    explicit ZmwChunkedFastqReaderPrivate(const std::string& fn, const size_t numChunks)
        : reader_{MakeFastqReaderImpl(std::move(fn), numChunks)}
    {
        assert(reader_->chunker_.NumChunks() != 0);
        Chunk(0);
    }

    void Chunk(size_t chunkId)
    {
        const auto& chunk = reader_->chunker_.Chunk(chunkId);
        remaining = chunk.NumRecords;
        reader_->Seek(chunk.FirstSeqOffset);
        currentChunkId_ = chunkId;
        firstRecord = true;
    }

    bool GetNext(FastqSequence& record)
    {
        if (remaining == 0) return false;
        record = reader_->ReadNextFastq(firstRecord);
        if (firstRecord) {
            record.Name(reader_->chunker_.Chunk(currentChunkId_).FirstSeqName);
            firstRecord = false;
        }
        --remaining;
        return true;
    }

    // reader
    std::unique_ptr<ZmwChunkedFastxReaderImpl> reader_;
    size_t currentChunkId_ = 0;
    bool firstRecord;
    size_t remaining;
};

static_assert(!std::is_copy_constructible<ZmwChunkedFastqReader>::value,
              "ZmwChunkedFastqReader(const ZmwChunkedFastqReader&) is not = delete");
static_assert(!std::is_copy_assignable<ZmwChunkedFastqReader>::value,
              "ZmwChunkedFastqReader& operator=(const ZmwChunkedFastqReader&) is not = delete");

ZmwChunkedFastqReader::ZmwChunkedFastqReader(const std::string& fn, const size_t numChunks)
    : internal::QueryBase<FastqSequence>{}
    , d_{std::make_unique<ZmwChunkedFastqReaderPrivate>(fn, numChunks)}
{
}

ZmwChunkedFastqReader::ZmwChunkedFastqReader(ZmwChunkedFastqReader&&) noexcept = default;

ZmwChunkedFastqReader& ZmwChunkedFastqReader::operator=(ZmwChunkedFastqReader&&) noexcept = default;

ZmwChunkedFastqReader::~ZmwChunkedFastqReader() = default;

size_t ZmwChunkedFastqReader::NumChunks() const { return d_->reader_->chunker_.NumChunks(); }

ZmwChunkedFastqReader& ZmwChunkedFastqReader::Chunk(size_t chunkId)
{
    d_->Chunk(chunkId);
    return *this;
}

size_t ZmwChunkedFastqReader::Chunk() const { return d_->currentChunkId_; }

bool ZmwChunkedFastqReader::GetNext(FastqSequence& record) { return d_->GetNext(record); }

}  // namespace BAM
}  // namespace PacBio
