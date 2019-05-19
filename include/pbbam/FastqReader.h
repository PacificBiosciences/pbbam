// File Description
/// \file FastqReader.h
/// \brief Defines the FastqReader class.
//
// Author: Derek Barnett

#ifndef FASTQREADER_H
#define FASTQREADER_H

#include <memory>
#include <string>
#include <vector>

#include "pbbam/FastqSequence.h"

#include "pbbam/internal/QueryBase.h"

namespace PacBio {
namespace BAM {

///
/// \brief The FastqReader provides sequential access to Fastq records.
///
class FastqReader : public internal::QueryBase<FastqSequence>
{
public:
    ///
    /// \brief Reads all Fastq sequences from a file
    ///
    /// \param fn   Fastq filename
    /// \return vector of FastqSequence results
    ///
    static std::vector<FastqSequence> ReadAll(const std::string& fn);

public:
    /// \name Constructors & Related Methods
    /// \{

    explicit FastqReader(const std::string& fn);

    FastqReader(const FastqReader&) = delete;
    FastqReader(FastqReader&&) noexcept;
    FastqReader& operator=(const FastqReader&) = delete;
    FastqReader& operator=(FastqReader&&) noexcept;
    ~FastqReader();

    /// \}

public:
    /// \name Sequence Access
    /// \{

    ///
    /// \brief GetNext
    ///
    /// \code{cpp}
    ///
    /// FastqReader reader{ fn };
    /// FastqSequence f;
    /// while (reader.GetNext(f)) {
    ///     // do stuff with f
    /// }
    /// \endcode
    ///
    /// \param[out] record
    /// \return success/failure
    ///
    bool GetNext(FastqSequence& record);

    /// \}

private:
    class FastqReaderPrivate;
    std::unique_ptr<FastqReaderPrivate> d_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // FASTQREADER_H
