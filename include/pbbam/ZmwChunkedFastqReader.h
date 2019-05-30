

// File Description
/// \file ZmwChunkedFastqReader.h
/// \brief Defines the ZmwChunkedFastqReader class.
//
// Author: Derek Barnett

#ifndef ZMWCHUNKEDFASTQREADER_H
#define ZMWCHUNKEDFASTQREADER_H

#include "pbbam/Config.h"

#include <memory>
#include <string>
#include <vector>

#include "pbbam/FastqSequence.h"
#include "pbbam/internal/QueryBase.h"

namespace PacBio {
namespace BAM {

///
/// \brief The ZmwChunkedFastqReader provides sequential access to FASTQ records,
///         where iteration is bounded by chunks of (unique) ZMWs.
///
class ZmwChunkedFastqReader : public internal::QueryBase<FastqSequence>
{
public:
    ///
    /// Construct a new ZMW-chunked FASTQ reader.
    ///
    /// \param fn           FASTQ file, must have a *.fai index
    /// \param numChunks    desired number of chunks
    ///
    /// Actual chunk count may be smaller than the requested number, if the input
    /// size is smaller.
    ///
    ZmwChunkedFastqReader(const std::string& fn, const size_t numChunks);

    ZmwChunkedFastqReader(ZmwChunkedFastqReader&&) noexcept;
    ZmwChunkedFastqReader& operator=(ZmwChunkedFastqReader&&) noexcept;
    ~ZmwChunkedFastqReader();

    ///
    /// \returns the number of chunks available.
    ///
    /// Actual chunk count may be smaller than the requested number, if the input
    /// size is smaller.
    ///
    size_t NumChunks() const;

    ///
    /// Sets current chunk to start iterating over.
    ///
    ZmwChunkedFastqReader& Chunk(size_t chunkId);

    ///
    /// \returns the current chunk in use
    ///
    size_t Chunk() const;

public:
    ///
    /// \brief GetNext
    ///
    /// Allows iteration with range-for:
    /// \code{cpp}
    ///
    /// ZmwChunkedFastqReader reader{fn, numChunks};
    /// reader.Chunk(4);
    /// for (const FastqSequence& seq : reader) {
    ///     // do stuff with seq
    /// }
    /// \endcode
    ///
    /// or you can iterate 'manually':
    /// \code{cpp}
    ///
    /// ZmwChunkedFastqReader reader{fn, numChunks};
    /// reader.Chunk(4);
    /// FastqSequence seq;
    /// while (reader.GetNext(seq)) {
    ///     // do stuff with seq
    /// }
    /// \endcode
    ///
    /// \param[out] record
    /// \return success/failure
    ///
    bool GetNext(FastqSequence& record);

private:
    class ZmwChunkedFastqReaderPrivate;
    std::unique_ptr<ZmwChunkedFastqReaderPrivate> d_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // ZMWCHUNKEDFASTQREADER_H
