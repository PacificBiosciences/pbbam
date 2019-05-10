// File Description
/// \file FastaReader.h
/// \brief Defines the FastaReader class.
//
// Author: Derek Barnett

#ifndef FASTAREADER_H
#define FASTAREADER_H

#include <memory>
#include <string>
#include <vector>

#include "pbbam/FastaSequence.h"

#include "pbbam/internal/QueryBase.h"

namespace PacBio {
namespace BAM {

///
/// \brief The FastaReader provides sequential access to FASTA records.
///
class FastaReader : public internal::QueryBase<FastaSequence>
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
    FastaReader(FastaReader&&);
    FastaReader& operator=(const FastaReader&) = delete;
    FastaReader& operator=(FastaReader&&);
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
    class FastaReaderPrivate;
    std::unique_ptr<FastaReaderPrivate> d_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // FASTAREADER_H
