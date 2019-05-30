// File Description
/// \file FastaReader.h
/// \brief Defines the FastaReader class.
//
// Author: Derek Barnett

#ifndef FASTAREADER_H
#define FASTAREADER_H

#include "pbbam/Config.h"

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

    FastaReader(FastaReader&&) noexcept;
    FastaReader& operator=(FastaReader&&) noexcept;
    ~FastaReader();

    /// \}

public:
    ///
    /// \brief GetNext
    ///
    /// Allows iteration with range-for:
    /// \code{cpp}
    ///
    /// FastaReader reader{fn};
    /// for (const FastaSequence& seq : reader) {
    ///     // do stuff with seq
    /// }
    /// \endcode
    ///
    /// or you can iterate 'manually':
    /// \code{cpp}
    ///
    /// FastaReader reader{fn};
    /// FastaSequence seq;
    /// while (reader.GetNext(seq)) {
    ///     // do stuff with seq
    /// }
    /// \endcode
    ///
    /// \param[out] record
    /// \return success/failure
    ///
    bool GetNext(FastaSequence& record);

private:
    class FastaReaderPrivate;
    std::unique_ptr<FastaReaderPrivate> d_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // FASTAREADER_H
