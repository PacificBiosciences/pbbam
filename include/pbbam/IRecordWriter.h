// File Description
/// \file IRecordWriter.h
/// \brief Defines the IRecordWriter interface.
//
// Author: Derek Barnett

#ifndef IRECORDWRITER_H
#define IRECORDWRITER_H

namespace PacBio {
namespace BAM {

class BamRecord;
class BamRecordImpl;

class IRecordWriter
{
public:
    virtual ~IRecordWriter() = default;

public:
    /// \brief Try to flush any buffered data to file.
    ///
    /// \note The underlying implementation may not necessarily flush buffered
    ///       data immediately, especially in a multithreaded writer situation.
    ///       Let the writer go out of scope to fully ensure flushing.
    ///
    /// \throws std::runtime_error if flush fails
    ///
    virtual void TryFlush() = 0;

    /// \brief Write a record to the output %BAM file.
    ///
    /// \param[in] record BamRecord object
    ///
    /// \throws std::runtime_error on failure to write
    ///
    virtual void Write(const BamRecord& record) = 0;

    /// \brief Write a record to the output %BAM file.
    ///
    /// \param[in] recordImpl BamRecordImpl object
    ///
    /// \throws std::runtime_error on failure to write
    ///
    virtual void Write(const BamRecordImpl& recordImpl) = 0;

protected:
    IRecordWriter() = default;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // IRECORDWRITER_H
