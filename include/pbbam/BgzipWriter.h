// File Description
/// \file BgzipWriter.h
/// \brief Defines the BgzipWriter class.
//
// Author: Derek Barnett

#ifndef BGZIPWRITER_H
#define BGZIPWRITER_H

#include "pbbam/Config.h"

#include <memory>
#include <string>

namespace PacBio {
namespace BAM {

struct BgzipWriterConfig
{
    // Compression level, equivalent to zlib-defined levels
    const size_t CompressionLevel = 0;

    /// Nnumber of threads for compression. If set to 0, the writer will attempt
    /// to determine a reasonable estimate. If set to 1, this will force
    /// single-threaded execution. No checks are made against an upper limit.
    const size_t NumThreads = 4;

    /// If true, write to <filename>.tmp, and rename to <filename> on closing.
    /// This provides for downstream checks to see if the file may be truncated
    /// due to early termination (e.g. a thrown exception).
    const bool UseTempFile = true;
};

/// \brief The BgzipWriter writes BGZF-compressed data to a file.
///
class PBBAM_EXPORT BgzipWriter
{
public:
    ///
    /// Create a BgzipWriter, using default configuration parameters.
    ///
    explicit BgzipWriter(std::string filename);

    ///
    /// Create a BgzipWriter, using configuration provided.
    ///
    BgzipWriter(std::string filename, const BgzipWriterConfig& config);

    BgzipWriter(BgzipWriter&&) noexcept;
    BgzipWriter& operator=(BgzipWriter&&) noexcept;
    ~BgzipWriter();

public:
    ///
    /// \brief Writes raw bytes to BGZF file.
    ///
    /// \param data         data buffer
    /// \param numBytes     num bytes (data size * sizeof(T))
    ///
    /// \returns number of bytes written
    ///
    size_t Write(const void* data, size_t numBytes);

    ///
    /// \brief Writes string data to BGZF file.
    ///
    /// \param data         data string
    ///
    /// \returns number of bytes written
    ///
    size_t Write(const std::string& data);

private:
    class BgzipWriterPrivate;
    std::unique_ptr<BgzipWriterPrivate> d_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // BGZIPWRITER_H
