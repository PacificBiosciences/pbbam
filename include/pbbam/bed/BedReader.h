#ifndef PBBAM_BED_BEDREADER_H
#define PBBAM_BED_BEDREADER_H

#include <pbbam/Config.h>

#include <pbbam/internal/QueryBase.h>

#include <pbcopper/data/GenomicInterval.h>

#include <memory>
#include <string>
#include <vector>

namespace PacBio {
namespace BED {

///
/// \brief The BedReader provides sequential access to BED records.
///
/// Supports plain text or gzipped (gzip or bgzip).
///
class BedReader : public BAM::internal::QueryBase<Data::GenomicInterval>
{
public:
    ///
    /// \brief Reads all BED intervals from a file
    ///
    /// \param fn   BED filename
    /// \return vector of intervals
    ///
    static std::vector<Data::GenomicInterval> ReadAll(const std::string& fn);

public:
    /// \name Constructors & Related Methods
    /// \{

    explicit BedReader(const std::string& fn);

    BedReader(BedReader&&) noexcept;
    BedReader& operator=(BedReader&&) noexcept;
    ~BedReader();

    /// \}

public:
    const std::string& Filename() const;

    ///
    /// \brief GetNext
    ///
    /// Allows iteration with range-for:
    /// \code{cpp}
    ///
    /// BedReader reader{fn};
    /// for (const auto& interval : reader) {
    ///     // do stuff with seq
    /// }
    /// \endcode
    ///
    /// or you can iterate 'manually':
    /// \code{cpp}
    ///
    /// BedReader reader{fn};
    /// GenomicInterval interval;
    /// while (reader.GetNext(interval)) {
    ///     // do stuff with interval
    /// }
    /// \endcode
    ///
    /// \param[out] record
    /// \return success/failure
    ///
    bool GetNext(Data::GenomicInterval& interval);

private:
    class BedReaderPrivate;
    std::unique_ptr<BedReaderPrivate> d_;
};

}  // namespace BED
}  // namespace PacBio

#endif  // PBBAM_BED_BEDREADER_H
