#ifndef PBBAM_INDEXEDBAMWRITER_H
#define PBBAM_INDEXEDBAMWRITER_H

#include <pbbam/Config.h>

#include <memory>
#include <string>

#include <pbbam/BamWriter.h>
#include <pbbam/IRecordWriter.h>
#include <pbbam/PbiBuilder.h>

namespace PacBio {
namespace BAM {

struct IndexedBamWriterConfig
{
    std::string outputFilename;
    BamHeader header;

    BamWriter::CompressionLevel bamCompressionLevel = BamWriter::DefaultCompression;
    PbiBuilder::CompressionLevel pbiCompressionLevel = PbiBuilder::DefaultCompression;

    // Number of threads used while writing to BAM file
    size_t numBamThreads = 4;
    // Number of threads used while writing to pbi file
    size_t numPbiThreads = 4;
    // Number of threads used while doing a trailing read of the BaM file being
    // written (to help compute indexes)
    size_t numGziThreads = 4;

    // Max size in memory for temporary files before flushing to disk.
    size_t tempFileBufferSize = 0x10000;
};
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
    /// \param[in] numGziThreads        number of threads used by the trailing
    ///                                 reader process used to help compute indexes.
    ///                                 If set to 0, the writer will attempt to
    ///                                 determine a reasonable estimate. If set
    ///                                 to 1, this will force single-threaded
    ///                                 execution. No checks are made against an
    ///                                 upper limit.
    /// \param[in] tempFileBufferBytes  Maximum number of bytes various temporary
    ///                                 files can use before they flush to disk.
    ///                                 Larger numbers require more resources but
    ///                                 may increase disk IO efficiency.
    ///
    /// \throws std::runtime_error if there was a problem
    ///
    IndexedBamWriter(
        const std::string& outputFilename, const BamHeader& header,
        const BamWriter::CompressionLevel bamCompressionLevel = BamWriter::DefaultCompression,
        const size_t numBamThreads = 4,
        const PbiBuilder::CompressionLevel pbiCompressionLevel = PbiBuilder::DefaultCompression,
        const size_t numPbiThreads = 4, const size_t numGziThreads = 4,
        const size_t tempFileBufferSize = 0x10000);

    /// \brief IndexedBamWRiter
    ///
    /// \param[in] config  Struct containing all the parameters used to construct
    ///                    this object.  See documentation for other constructor
    ///                    for more details
    IndexedBamWriter(const IndexedBamWriterConfig& config)
        : IndexedBamWriter(config.outputFilename, config.header, config.bamCompressionLevel,
                           config.numBamThreads, config.pbiCompressionLevel, config.numPbiThreads,
                           config.numGziThreads, config.tempFileBufferSize)
    {
    }

    IndexedBamWriter(IndexedBamWriter&&) noexcept;
    IndexedBamWriter& operator=(IndexedBamWriter&&) noexcept;
    ~IndexedBamWriter();

public:
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
    class IndexedBamWriterPrivate2;
    std::unique_ptr<IndexedBamWriterPrivate2> d_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_INDEXEDBAMWRITER_H
