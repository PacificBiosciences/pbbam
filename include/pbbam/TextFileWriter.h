
// File Description
/// \file TextFileWriter.h
/// \brief Defines the TextFileWriter class.
//
// Author: Derek Barnett

#ifndef PBBAM_TEXTFILEWRITER_H
#define PBBAM_TEXTFILEWRITER_H

#include "pbbam/Config.h"

#include <memory>
#include <string>

namespace PacBio {
namespace BAM {

///
/// \brief The TextFileWriter class provides line-by-line writing text files.
///
/// Supports plain text or gzipped. For explicitly-bgzipped text, use
/// BgzipWriter instead.
///
/// \note This is a general-purpose file writer. For FASTA/FASTQ, use the
///       dedicated FastaReader/FastqReader or BgzipFastaWriter/BgzipFastaWriter
///       for better performance.
///
class TextFileWriter
{
public:
    ///
    /// \brief TextLineReader
    ///
    /// \param filename  suffix ".gz" indicates gzipped output
    ///
    explicit TextFileWriter(const std::string& filename);

    TextFileWriter(TextFileWriter&&) noexcept;
    TextFileWriter& operator=(TextFileWriter&&) noexcept;
    ~TextFileWriter();

public:
    ///
    /// \brief Write
    ///
    ///
    /// \param line
    ///
    void Write(const std::string& line);

private:
    class TextFileWriterPrivate;
    std::unique_ptr<TextFileWriterPrivate> d_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_TEXTFILEWRITER_H