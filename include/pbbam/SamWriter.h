// File Description
/// \file SamWriter.h
/// \brief Defines the SamWriter class.
//
// Author: Derek Barnett

#ifndef SAMWRITER_H
#define SAMWRITER_H

#include <memory>
#include <string>
#include "pbbam/BamHeader.h"
#include "pbbam/BamRecord.h"
#include "pbbam/IRecordWriter.h"

namespace PacBio {
namespace BAM {

namespace internal {
class SamWriterPrivate;
}

/// \brief The SamWriter class provides a writing interface for creating
///        new SAM files.
///
/// \note The underlying buffered data may not be flushed to the file until the
///       destructor is called. Trying to access the file (reading, stat-ing,
///       indexing, etc.) before the SamWriter is destroyed yields undefined
///       behavior. Enclose the SamWriter in some form of local scope (curly
///       braces, a separate function, etc.) to ensure that its destructor is
///       called before proceeding to read-based operations.
///
/// \code{.cpp}
///  {
///     SamWriter w(...);
///     // write data
///  }
///  // now safe to access the new file
/// \endcode
///
///
class SamWriter : public IRecordWriter
{
public:
    /// \brief Opens a SAM file for writing & writes the header information.
    ///
    /// \note Set \p filename to "-" for stdout.
    ///
    /// \param[in] filename     path to output SAM file
    /// \param[in] header       BamHeader object
    ///
    /// \throws std::runtime_error if there was a problem opening the file for
    ///         writing or if an error occurred while writing the header
    ///
    SamWriter(std::string filename, const BamHeader& header);

    /// Fully flushes all buffered data & closes file.
    ///
    ~SamWriter() override;

    SamWriter(const SamWriter&) = delete;
    SamWriter(SamWriter&&) = default;
    SamWriter& operator=(const SamWriter&) = delete;
    SamWriter& operator=(SamWriter&&) = default;

public:
    /// \brief Try to flush any buffered data to file.
    ///
    /// \note The underlying implementation may not necessarily flush buffered
    ///       data immediately, especially in a multithreaded writer situation.
    ///       Let the SamWriter go out of scope to fully ensure flushing.
    ///
    /// \throws std::runtime_error if flush fails
    ///
    void TryFlush() override;

    /// \brief Write a record to the output SAM file.
    ///
    /// \param[in] record BamRecord object
    ///
    /// \throws std::runtime_error on failure to write
    ///
    void Write(const BamRecord& record) override;

    /// \brief Write a record to the output SAM file.
    ///
    /// \param[in] recordImpl BamRecordImpl object
    ///
    /// \throws std::runtime_error on failure to write
    ///
    void Write(const BamRecordImpl& recordImpl) override;

private:
    std::unique_ptr<internal::SamWriterPrivate> d_;
};

}  // namesapce BAM
}  // namespace PacBio

#endif  // SAMWRITER_H
