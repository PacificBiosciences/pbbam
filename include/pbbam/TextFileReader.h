#ifndef PBBAM_TEXTFILEREADER_H
#define PBBAM_TEXTFILEREADER_H

#include <pbbam/Config.h>

#include <memory>
#include <string>

#include <pbbam/internal/QueryBase.h>

namespace PacBio {
namespace BAM {

///
/// \brief The TextFileReader class provides line-by-line access to text files.
///
/// Supports plain text or gzipped (gzip or bgzip).
///
/// \note This is a general-purpose file reader. For FASTA/FASTQ, use the dedicated
///       FastaReader or FastqReader for better performance.
///
class TextFileReader : public internal::QueryBase<std::string>
{
public:
    ///
    /// \brief Reads all lines from a text file
    ///
    /// \param fn    filename
    /// \return vector of lines
    ///
    static std::vector<std::string> ReadAll(const std::string& fn);

public:
    ///
    /// \brief TextLineReader
    ///
    /// \param filename
    ///
    explicit TextFileReader(std::string filename);

    TextFileReader(TextFileReader&&) noexcept;
    TextFileReader& operator=(TextFileReader&&) noexcept;
    ~TextFileReader();

public:
    const std::string& Filename() const;

    //
    /// \brief GetNext
    ///
    /// Allows iteration with range-for:
    /// \code{cpp}
    ///
    /// TextFileReader reader{fn};
    /// for (const std::string& line : reader) {
    ///     // do stuff with line
    /// }
    /// \endcode
    ///
    /// or you can iterate 'manually':
    /// \code{cpp}
    ///
    /// TextFileReader reader{fn};
    /// std::string line;
    /// while (reader.GetNext(line)) {
    ///     // do stuff with line
    /// }
    /// \endcode
    ///
    /// \param[out] line
    /// \return success/failure
    ///
    bool GetNext(std::string& line) override;

private:
    class TextFileReaderPrivate;
    std::unique_ptr<TextFileReaderPrivate> d_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_TEXTFILEREADER_H
