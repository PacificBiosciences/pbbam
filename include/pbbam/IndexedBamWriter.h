// File Description
/// \file IndexedBamWriter.h
/// \brief Defines the IndexedBamWriter class.
//
// Author: Derek Barnett

#ifndef INDEXEDBAMWRITER_H
#define INDEXEDBAMWRITER_H

#include <memory>
#include <string>

#include "pbbam/BamWriter.h"
#include "pbbam/Config.h"
#include "pbbam/IRecordWriter.h"
#include "pbbam/PbiBuilder.h"

namespace PacBio {
namespace BAM {

class BamHeader;
class BamRecord;
class BamRecordImpl;

namespace internal {
class IndexedBamWriterPrivate2;
}

///
/// \brief The IndexedBamWriter class
///
///
///
///
///
///
class IndexedBamWriter : public IRecordWriter
{
public:
    ///
    /// \brief IndexedBamWriter
    ///
    /// \param[in] filename             path to output %BAM file
    /// \param[in] header               BAM file header
    ///
    /// \param[in] bamCompressionLevel  zlib compression level for output BAM
    /// \param[in] numBamThreads        number of threads for BAM compression.
    ///                                 If set to 0, the writer will attempt to
    ///                                 determine a reasonable estimate. If set
    ///                                 to 1, this will force single-threaded
    ///                                 execution. No checks are made against an
    ///                                 upper limit.
    ///
    /// \param[in] pbiCompressionLevel  zlib compression level for output PBI
    /// \param[in] numPbiThreads        number of threads for PBI compression.
    ///                                 If set to 0, the writer will attempt to
    ///                                 determine a reasonable estimate. If set
    ///                                 to 1, this will force single-threaded
    ///                                 execution. No checks are made against an
    ///                                 upper limit.
    ///
    /// \throws std::runtime_error if there was a problem
    ///
    IndexedBamWriter(
        const std::string& outputFilename, const BamHeader& header,
        const BamWriter::CompressionLevel bamCompressionLevel = BamWriter::DefaultCompression,
        const size_t numBamThreads = 4,
        const PbiBuilder::CompressionLevel pbiCompressionLevel = PbiBuilder::DefaultCompression,
        const size_t numPbiThreads = 4);

    IndexedBamWriter(const IndexedBamWriter&) = delete;
    IndexedBamWriter(IndexedBamWriter&&) = default;
    IndexedBamWriter& operator=(const IndexedBamWriter&) = delete;
    IndexedBamWriter& operator=(IndexedBamWriter&&) = default;
    ~IndexedBamWriter() override;

public:
    ///
    /// \brief TryFlush
    ///
    void TryFlush() override;

    ///
    /// \brief Write
    ///
    /// \param[in] record
    ///
    void Write(const BamRecord& record) override;

    ///
    /// \brief Write
    ///
    /// \param[in] record
    ///
    void Write(const BamRecordImpl& record) override;

private:
    std::unique_ptr<internal::IndexedBamWriterPrivate2> d_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // INDEXEDBAMWRITER_H
