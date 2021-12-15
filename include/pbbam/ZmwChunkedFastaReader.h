#ifndef PBBAM_ZMWCHUNKEDFASTAREADER_H
#define PBBAM_ZMWCHUNKEDFASTAREADER_H

#include <pbbam/Config.h>

#include <pbbam/FastaSequence.h>
#include <pbbam/internal/QueryBase.h>

#include <memory>
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

///
/// \brief The ZmwChunkedFastaReader provides sequential access to FASTA records,
///         where iteration is bounded by chunks of (unique) ZMWs.
///
class ZmwChunkedFastaReader : public internal::QueryBase<FastaSequence>
{
public:
    ///
    /// Construct a new ZMW-chunked FASTA reader.
    ///
    /// \param fn           FASTA file, must have a *.fai index
    /// \param numChunks    desired number of chunks
    ///
    /// Actual chunk count may be smaller than the requested number, if the input
    /// size is smaller.
    ///
    ZmwChunkedFastaReader(const std::string& fn, size_t numChunks);

    ZmwChunkedFastaReader(ZmwChunkedFastaReader&&) noexcept;
    ZmwChunkedFastaReader& operator=(ZmwChunkedFastaReader&&) noexcept;
    ~ZmwChunkedFastaReader();

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
    ZmwChunkedFastaReader& Chunk(size_t chunkId);

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
    /// ZmwChunkedFastaReader reader{fn, numChunks};
    /// reader.Chunk(4);
    /// for (const FastaSequence& seq : reader) {
    ///     // do stuff with seq
    /// }
    /// \endcode
    ///
    /// or you can iterate 'manually':
    /// \code{cpp}
    ///
    /// ZmwChunkedFastaReader reader{fn, numChunks};
    /// reader.Chunk(4);
    /// FastaSequence seq;
    /// while (reader.GetNext(seq)) {
    ///     // do stuff with seq
    /// }
    /// \endcode
    ///
    /// \param[out] record
    /// \return success/failure
    ///
    bool GetNext(FastaSequence& record);

private:
    class ZmwChunkedFastaReaderPrivate;
    std::unique_ptr<ZmwChunkedFastaReaderPrivate> d_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_ZMWCHUNKEDFASTAREADER_H
