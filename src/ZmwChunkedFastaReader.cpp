// File Description
/// \file ZmwChunkedFastaReader.cpp
/// \brief Implements the ZmwChunkedFastaReader class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/ZmwChunkedFastaReader.h"

#include <cstdio>

#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include <htslib/kseq.h>

#include "pbbam/FaiIndex.h"
#include "pbbam/FormatUtils.h"

#include "MemoryUtils.h"
#include "ZmwChunkedFastxBgzfReader.h"
#include "ZmwChunkedFastxReaderImpl.h"
#include "ZmwChunkedFastxTextReader.h"

namespace PacBio {
namespace BAM {
namespace {

std::unique_ptr<ZmwChunkedFastxReaderImpl> MakeFastaReaderImpl(std::string filename,
                                                               const size_t numChunks)
{
    // validate extension
    if (!FormatUtils::IsFastaFilename(filename)) {
        throw std::runtime_error{"ZmwChunkedFastaReader: filename '" + filename +
                                 "' is not recognized as a FASTA file."};
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
            msg << "ZmwChunkedFastaReader: random-access is not supported for plain gzipped "
                   "file "
                << filename << "\n\n"
                << "Compressed files must be bgzipped, with accompanying *.gzi "
                   "index.\n\n"
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
            throw std::runtime_error{
                "ZmwChunkedFastaReader: could not determine compression type for file: " +
                filename};
    }
}

}  // namespace

class ZmwChunkedFastaReader::ZmwChunkedFastaReaderPrivate
{
public:
    explicit ZmwChunkedFastaReaderPrivate(const std::string& fn, const size_t numChunks)
        : reader_{MakeFastaReaderImpl(std::move(fn), numChunks)}
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

    bool GetNext(FastaSequence& record)
    {
        if (remaining == 0) return false;
        record = reader_->ReadNextFasta(firstRecord);
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

static_assert(!std::is_copy_constructible<ZmwChunkedFastaReader>::value,
              "ZmwChunkedFastaReader(const ZmwChunkedFastaReader&) is not = delete");
static_assert(!std::is_copy_assignable<ZmwChunkedFastaReader>::value,
              "ZmwChunkedFastaReader& operator=(const ZmwChunkedFastaReader&) is not = delete");

ZmwChunkedFastaReader::ZmwChunkedFastaReader(const std::string& fn, const size_t numChunks)
    : internal::QueryBase<FastaSequence>{}
    , d_{std::make_unique<ZmwChunkedFastaReaderPrivate>(fn, numChunks)}
{
}

ZmwChunkedFastaReader::ZmwChunkedFastaReader(ZmwChunkedFastaReader&&) noexcept = default;

ZmwChunkedFastaReader& ZmwChunkedFastaReader::operator=(ZmwChunkedFastaReader&&) noexcept = default;

ZmwChunkedFastaReader::~ZmwChunkedFastaReader() = default;

size_t ZmwChunkedFastaReader::NumChunks() const { return d_->reader_->chunker_.NumChunks(); }

ZmwChunkedFastaReader& ZmwChunkedFastaReader::Chunk(size_t chunkId)
{
    d_->Chunk(chunkId);
    return *this;
}

size_t ZmwChunkedFastaReader::Chunk() const { return d_->currentChunkId_; }

bool ZmwChunkedFastaReader::GetNext(FastaSequence& record) { return d_->GetNext(record); }

}  // namespace BAM
}  // namespace PacBio
