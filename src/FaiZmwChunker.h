// File Description
/// \file FaiZmwChunker.h
/// \brief Defines the FaiZmwChunker enum.
//
// Author: Derek Barnett

#ifndef FAIZMWCHUNKER_H
#define FAIZMWCHUNKER_H

#include "pbbam/Config.h"

#include <string>
#include <vector>

#include "pbbam/FaiIndex.h"

namespace PacBio {
namespace BAM {

struct FaiZmwChunk
{
    /// Name of first entry
    std::string FirstSeqName;

    /// File offset to the sequence of the chunk's first entry.
    uint64_t FirstSeqOffset;

    // Total number of records in chunk.
    size_t NumRecords;

    // Number of unique ZMWs
    size_t NumZmws;
};

///
/// \brief The FaiZmwChunker takes a FAI index and bins unique ZMW hole numbers
///        into chunks.
///
class FaiZmwChunker
{
public:
    ///
    /// \brief Construct a new FaiZmwChunker
    ///
    /// \param index        FAI index
    /// \param numChunks    desired number of chunks
    ///
    /// Actual chunk count may be smaller than the requested number, if the input
    /// size is smaller.
    ///
    FaiZmwChunker(const FaiIndex& index, const size_t numChunks);

    ///
    /// \brief Construct a new FaiZmwChunker
    ///
    /// \param filename     FAI filename
    /// \param numChunks    desired number of chunks
    ///
    /// Actual chunk count may be smaller than the requested number, if the input
    /// size is smaller.
    ///
    FaiZmwChunker(const std::string& filename, const size_t numChunks);

    FaiZmwChunker(const FaiZmwChunker&);
    FaiZmwChunker(FaiZmwChunker&&) noexcept;
    FaiZmwChunker& operator=(const FaiZmwChunker&);
    FaiZmwChunker& operator=(FaiZmwChunker&&) noexcept;
    ~FaiZmwChunker();

public:
    const FaiZmwChunk& Chunk(size_t chunk) const;

    size_t NumChunks() const;

    // size_t MaxChunkSize(size_t index) const;

private:
    std::vector<FaiZmwChunk> chunks_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // FAIZMWCHUNKER_H
