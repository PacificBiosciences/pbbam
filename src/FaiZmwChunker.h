#ifndef PBBAM_FAIZMWCHUNKER_H
#define PBBAM_FAIZMWCHUNKER_H

#include <pbbam/Config.h>

#include <pbbam/FaiIndex.h>

#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

struct FaiZmwChunk
{
    /// Name of first entry
    std::string FirstSeqName;

    /// File offset to the sequence of the chunk's first entry.
    uint64_t FirstSeqOffset;

    // Total number of records in chunk.
    std::size_t NumRecords;

    // Number of unique ZMWs
    std::size_t NumZmws;
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
    FaiZmwChunker(const FaiIndex& index, std::size_t numChunks);

    ///
    /// \brief Construct a new FaiZmwChunker
    ///
    /// \param filename     FAI filename
    /// \param numChunks    desired number of chunks
    ///
    /// Actual chunk count may be smaller than the requested number, if the input
    /// size is smaller.
    ///
    FaiZmwChunker(const std::string& filename, std::size_t numChunks);

    FaiZmwChunker(const FaiZmwChunker&);
    FaiZmwChunker(FaiZmwChunker&&) noexcept;
    FaiZmwChunker& operator=(const FaiZmwChunker&);
    FaiZmwChunker& operator=(FaiZmwChunker&&) noexcept;
    ~FaiZmwChunker();

public:
    const FaiZmwChunk& Chunk(std::size_t chunk) const;
    std::size_t NumChunks() const;

private:
    std::vector<FaiZmwChunk> chunks_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_FAIZMWCHUNKER_H
