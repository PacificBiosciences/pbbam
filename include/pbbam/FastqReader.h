#ifndef PBBAM_FASTQREADER_H
#define PBBAM_FASTQREADER_H

#include <pbbam/Config.h>

#include <pbbam/FastqSequence.h>
#include <pbbam/internal/QueryBase.h>

#include <memory>
#include <string>
#include <vector>

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

    FastqReader(FastqReader&&) noexcept;
    FastqReader& operator=(FastqReader&&) noexcept;
    virtual ~FastqReader();

    /// \}

public:
    ///
    /// \brief GetNext
    ///
    /// Allows iteration with range-for:
    /// \code{cpp}
    ///
    /// FastqReader reader{fn};
    /// for (const FastqSequence& seq : reader) {
    ///     // do stuff with seq
    /// }
    /// \endcode
    ///
    /// or you can iterate 'manually':
    /// \code{cpp}
    ///
    /// FastqReader reader{fn};
    /// FastqSequence seq;
    /// while (reader.GetNext(seq)) {
    ///     // do stuff with seq
    /// }
    /// \endcode
    ///
    /// \param[out] record
    /// \return success/failure
    ///
    bool GetNext(FastqSequence& record);

private:
    class FastqReaderPrivate;
    std::unique_ptr<FastqReaderPrivate> d_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_FASTQREADER_H
