#ifndef PBBAM_PBIFILTERQUERY_H
#define PBBAM_PBIFILTERQUERY_H

#include <pbbam/Config.h>

#include <pbbam/PbiFilter.h>
#include <pbbam/internal/QueryBase.h>

#include <vector>

namespace PacBio {
namespace BAM {

/// \brief The PbiFilter class provides iterable access to a DataSet's %BAM
///        records, limiting results to those matching filter criteria.
///
/// Example:
/// \include code/PbiFilterQuery.txt
///
/// \note Currently, all %BAM files must have a corresponding ".pbi" index file.
///       Use BamFile::EnsurePacBioIndexExists before creating the query if one
///       may not be present.
///
class PBBAM_EXPORT PbiFilterQuery : public internal::IQuery
{
public:
    ///
    /// \brief Creates a new PbiFilterQuery, limiting record results to only
    ///        those matching filter criteria defined in the DataSet XML.
    ///
    /// \param[in] dataset input data source(s)
    ///
    /// \throws std::runtime_error on failure to open/read underlying %BAM or
    ///         PBI files.
    ///
    PbiFilterQuery(const DataSet& dataset);

    PbiFilterQuery(const DataSet& dataset, const PbiIndexCache& cache);

    /// \brief Creates a new PbiFilterQuery, limiting record results to only
    ///        those matching filter criteria
    ///
    /// \param[in] filter   filtering criteria
    /// \param[in] dataset  input data source(s)
    ///
    /// \throws std::runtime_error on failure to open/read underlying %BAM or
    ///         PBI files.
    ///
    PbiFilterQuery(const PbiFilter& filter, const DataSet& dataset);

    PbiFilterQuery(const PbiFilter& filter, const DataSet& dataset, const PbiIndexCache& cache);

    PbiFilterQuery(PbiFilterQuery&&) noexcept;
    PbiFilterQuery& operator=(PbiFilterQuery&&) noexcept;
    ~PbiFilterQuery() override;

    /// \brief Main iteration point for record access.
    ///
    /// Most client code should not need to use this method directly. Use
    /// iterators instead.
    ///
    bool GetNext(BamRecord& r) override;

    /// \brief Return number of records that pass the provided filter
    ///
    std::uint32_t NumReads() const;

private:
    class PbiFilterQueryPrivate;
    std::unique_ptr<PbiFilterQueryPrivate> d_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_PBIFILTERQUERY_H
