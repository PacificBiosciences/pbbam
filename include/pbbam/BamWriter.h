#ifndef PBBAM_BAMWRITER_H
#define PBBAM_BAMWRITER_H

#include <pbbam/Config.h>

#include <pbbam/BamFile.h>
#include <pbbam/BamHeader.h>
#include <pbbam/BamRecord.h>
#include <pbbam/IRecordWriter.h>

#include <htslib/sam.h>

#include <string>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace BAM {

/// \brief The BamWriter class provides a writing interface for creating
///        new %BAM files.
///
/// \note The underlying buffered data may not be flushed to the file until the
///       destructor is called. Trying to access the file (reading, stat-ing,
///       indexing, etc.) before the BamWriter is destroyed yields undefined
///       behavior. Enclose the BamWriter in some form of local scope (curly
///       braces, a separate function, etc.) to ensure that its destructor is
///       called before proceeding to read-based operations.
///
/// \code{.cpp}
///  {
///     BamWriter w(...);
///     // write data
///  }
///  // now safe to access the new file
/// \endcode
///
///
class PBBAM_EXPORT BamWriter : public IRecordWriter
{
public:
    /// \brief This enum allows you to control the compression level of the
    ///        output %BAM file.
    ///
    /// Values are equivalent to zlib compression levels. See its documentation
    /// for more details: http://www.zlib.net/manual.html
    ///
    enum CompressionLevel
    {
        CompressionLevel_0 = 0,
        CompressionLevel_1 = 1,
        CompressionLevel_2 = 2,
        CompressionLevel_3 = 3,
        CompressionLevel_4 = 4,
        CompressionLevel_5 = 5,
        CompressionLevel_6 = 6,
        CompressionLevel_7 = 7,
        CompressionLevel_8 = 8,
        CompressionLevel_9 = 9,

        DefaultCompression = -1,
        NoCompression = CompressionLevel_0,
        FastCompression = CompressionLevel_1,
        BestCompression = CompressionLevel_9
    };

    /// \brief This enum allows you to control whether BAI bin numbers are
    ///        calculated for output records.
    ///
    /// For most cases, the default behavior (ON) should be retained for maximum
    /// compatibility with downstream tools (e.g. samtools index). Disabling bin
    /// calculation should only be used if all records are known to never be
    /// mapped, and even then only if profiling revelas the calculation to
    /// affect extremely performance-sensitive, "critical paths".
    ///
    enum BinCalculationMode
    {
        BinCalculation_ON = 0,
        BinCalculation_OFF
    };

    ///
    /// \brief The Config struct provides a "parameter object" for BamWriter
    ///        settings. This allows for writer configuration without having to
    ///        refer to ordering of parameters, default values, etc.
    ///
    struct Config
    {
        // zlib compression level
        CompressionLevel compressionLevel = DefaultCompression;

        // The number of threads for compression. If set to 0, BamWriter will
        // attempt to determine a reasonable estimate. If set to 1, this will
        // force single-threaded execution. No checks are made against an upper limit.
        std::size_t numThreads = 4;

        // If ON, ensures that proper BAI bin numbers are provided for all records.
        BamWriter::BinCalculationMode binCalculationMode = BamWriter::BinCalculation_ON;

        // If true, write to <filename>.tmp, and rename  to <filename> in dtor.
        // This allows downstream checks to see if BAM file may be truncated
        // due to early termination (e.g. a thrown exception). If false, write
        // directly to <filename>.
        bool useTempFile = true;
    };

public:
    /// \name Constructors & Related Methods
    /// \{

    /// \brief Opens a %BAM file for writing & writes the header information.
    ///
    /// \note Set \p filename to "-" for stdout.
    ///
    /// \param[in] filename         path to output %BAM file
    /// \param[in] header           BamHeader object
    /// \param[in] compressionLevel zlib compression level
    /// \param[in] numThreads       number of threads for compression. If set to
    ///                             0, BamWriter will attempt to determine a
    ///                             reasonable estimate. If set to 1, this will
    ///                             force single-threaded execution. No checks
    ///                             are made against an upper limit.
    ///
    /// \param[in] binCalculationMode BAI bin calculation mode. The default
    ///            behavior will ensure proper bin numbers are provided for all
    ///            records written. This extra step may turned off when bin
    ///            numbers are not needed. Though if in doubt, keep the default.
    ///
    /// \param[in] useTempFile      If true, write to <filename>.tmp, and rename
    ///                             to <filename>. This provides for downstream
    ///                             checks to see if BAM file may be truncated
    ///                             due to early termination (a thrown exception).
    ///
    /// \throws std::runtmie_error if there was a problem opening the file for
    ///         writing or if an error occurred while writing the header
    ///
    BamWriter(const std::string& filename, const BamHeader& header,
              BamWriter::CompressionLevel compressionLevel = BamWriter::DefaultCompression,
              std::size_t numThreads = 4,
              BinCalculationMode binCalculationMode = BamWriter::BinCalculation_ON,
              bool useTempFile = true);

    ///
    /// \brief Opens a %BAM file for writing & writes the header information.
    ///
    /// \param[in] filename     path to output %BAM file
    /// \param[in] header       BamHeader object
    /// \param[in] config       container for add'l configuration options
    ///
    /// \throws std::runtmie_error if there was a problem opening the file for
    ///         writing or if an error occurred while writing the header
    ///
    BamWriter(const std::string& filename, const BamHeader& header,
              const BamWriter::Config& config);

    BamWriter(BamWriter&&) noexcept;
    BamWriter& operator=(BamWriter&&) noexcept;

    /// Fully flushes all buffered data & closes file.
    ~BamWriter() override;

    /// \}

public:
    /// \name Data Writing & Resource Management
    /// \{

    /// \brief Try to flush any buffered data to file.
    ///
    /// \note The underlying implementation doesn't necessarily flush buffered
    ///       data immediately, especially in a multithreaded writer situation.
    ///       Let the BamWriter go out of scope to fully ensure flushing.
    ///
    /// \throws std::runtime_error if flush fails
    ///
    void TryFlush() override;

    /// \brief Write a record to the output %BAM file.
    ///
    /// \param[in] record BamRecord object
    ///
    /// \throws std::runtime_error on failure to write
    ///
    void Write(const BamRecord& record) override;

    /// \brief Write a record to the output %BAM file.
    ///
    /// \param[in] record BamRecord object
    /// \param[out] vOffset BGZF virtual offset to start of \p record
    ///
    /// \throws std::runtime_error on failure to write
    ///
    void Write(const BamRecord& record, int64_t* vOffset);

    /// \brief Write a record to the output %BAM file.
    ///
    /// \param[in] recordImpl BamRecordImpl object
    ///
    /// \throws std::runtime_error on failure to write
    ///
    void Write(const BamRecordImpl& recordImpl) override;

    /// \}

private:
    class BamWriterPrivate;
    std::unique_ptr<BamWriterPrivate> d_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_BAMWRITER_H
