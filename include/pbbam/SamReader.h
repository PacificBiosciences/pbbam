// Author: Derek Barnett

#ifndef PBBAM_SAMREADER_H
#define PBBAM_SAMREADER_H

#include "pbbam/Config.h"

#include <memory>
#include <string>

#include "pbbam/BamHeader.h"
#include "pbbam/BamRecord.h"
#include "pbbam/internal/QueryBase.h"

namespace PacBio {
namespace BAM {

/// \brief The SamReader class provides basic read-access to a SAM file.
///
class PBBAM_EXPORT SamReader : public internal::IQuery
{
public:
    /// \brief Opens SAM for streaming from stdin
    SamReader();

    /// \brief Opens SAM file for reading.
    ///
    /// \param[in] fn   SAM filename
    /// \throw std::runtime_error if failed to open
    ///
    explicit SamReader(std::string fn);

    virtual ~SamReader();

public:
    /// \returns SAM filename
    const std::string& Filename() const;

    /// \returns BamHeader object from SAM
    const BamHeader& Header() const;

public:
    /// \brief Fetches the next record.
    ///
    /// \param[out] record  next BamRecord object. Should not be used if method
    ///                     returns false.
    ///
    /// \returns true if record was read successfully. Returns false if EOF (or
    ///          end of iterator in derived readers). False is not an error,
    ///          it indicates "end of data".
    ///
    /// \throw std::runtime_error if failed to read from file (e.g. possible
    ///        truncated or corrupted file).
    ///
    bool GetNext(BamRecord& record);

private:
    class SamReaderPrivate;
    std::unique_ptr<SamReaderPrivate> d_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_SAMREADER_H
