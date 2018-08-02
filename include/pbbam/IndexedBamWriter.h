// File Description
/// \file IndexedBamWriter.h
/// \brief Defines the IndexedBamWriter class.
//
// Author: Derek Barnett

#ifndef INDEXEDBAMWRITER_H
#define INDEXEDBAMWRITER_H

#include <memory>
#include <string>

#include "pbbam/Config.h"
#include "pbbam/IRecordWriter.h"

namespace PacBio {
namespace BAM {

class BamHeader;
class BamRecord;
class BamRecordImpl;

namespace internal {
class IndexedBamWriterPrivate;
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
    /// \param[in] filename         path to output %BAM file
    /// \param[in] header           BamHeader object
    ///
    /// \throws std::runtime_error if there was a problem opening the file for
    ///         writing or if an error occurred while writing the header
    ///
    IndexedBamWriter(const std::string& outputFilename, const BamHeader& header);

    ~IndexedBamWriter() override;

    IndexedBamWriter(const IndexedBamWriter&) = delete;
    IndexedBamWriter(IndexedBamWriter&&) = delete;
    IndexedBamWriter& operator=(const IndexedBamWriter&) = delete;
    IndexedBamWriter& operator=(IndexedBamWriter&&) = delete;

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
    std::unique_ptr<internal::IndexedBamWriterPrivate> d_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // INDEXEDBAMWRITER_H
