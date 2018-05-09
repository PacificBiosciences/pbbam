// File Description
/// \file FastaReader.h
/// \brief Defines the FastaReader class.
//
// Author: Derek Barnett

#ifndef FASTAREADER_H
#define FASTAREADER_H

#include <memory>
#include <vector>
#include "pbbam/FastaSequence.h"

namespace PacBio {
namespace BAM {

namespace internal {
struct FastaReaderPrivate;
}

///
/// \brief The FastaReader provides sequential access to FASTA records.
///
class FastaReader
{
public:
    ///
    /// \brief Reads all FASTA sequences from a file
    ///
    /// \param fn   FASTA filename
    /// \return vector of FastaSequence results
    ///
    static std::vector<FastaSequence> ReadAll(const std::string& fn);

public:
    /// \name Constructors & Related Methods
    /// \{

    explicit FastaReader(const std::string& fn);
    FastaReader(const FastaReader&) = delete;
    FastaReader(FastaReader&&) = default;
    FastaReader& operator=(const FastaReader&) = delete;
    FastaReader& operator=(FastaReader&&) = default;
    ~FastaReader();

    /// \}

public:
    /// \name Sequence Access
    /// \{

    ///
    /// \brief GetNext
    ///
    /// \code{cpp}
    ///
    /// FastaReader reader{ fn };
    /// FastaSequence f;
    /// while (reader.GetNext(f)) {
    ///     // do stuff with f
    /// }
    /// \endcode
    ///
    /// \param[out] record
    /// \return success/failure
    ///
    bool GetNext(FastaSequence& record);

    /// \}

private:
    std::unique_ptr<internal::FastaReaderPrivate> d_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // FASTAREADER_H
