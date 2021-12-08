#ifndef PBBAM_GENOMICINTERVALQUERY_H
#define PBBAM_GENOMICINTERVALQUERY_H

#include <pbbam/Config.h>

#include <pbbam/BaiIndexCache.h>
#include <pbbam/BamRecord.h>
#include <pbbam/DataSet.h>
#include <pbbam/internal/QueryBase.h>

#include <pbcopper/data/GenomicInterval.h>

#include <memory>

namespace PacBio {
namespace BAM {

/// \brief The GenomicIntervalQuery class provides iterable access to a
///        DataSet's %BAM records, limiting results to those overlapping a
///        GenomicInterval.
///
/// Example:
/// \include code/GenomicIntervalQuery.txt
///
/// \note Currently, all %BAM files must have a corresponding ".bai" index file.
///       Use BamFile::EnsureStandardIndexExists before creating the query if
///       one may not be present.
///
class PBBAM_EXPORT GenomicIntervalQuery : public internal::IQuery
{
public:
    /// \brief Constructs a new GenomiIntervalQuery, that can be used to retrieve
    ///        only those records overlapping a GenomicInterval.
    ///
    /// \note Using this constructor means that an interval must be provided, via
    ///       query.Interval(i), before iterating.
    ///
    /// \param[in] dataset  input data source(s)
    ///
    /// \throws std::runtime_error on failure to open/read underlying %BAM or
    ///         BAI files.
    ///
    GenomicIntervalQuery(const DataSet& dataset);
    GenomicIntervalQuery(const DataSet& dataset, const BaiIndexCache& cache);

    /// \brief Constructs a new GenomiIntervalQuery, limiting record results to
    ///        only those overalpping a GenomicInterval.
    ///
    /// \param[in] interval genomic interval of interest
    /// \param[in] dataset  input data source(s)
    ///
    /// \throws std::runtime_error on failure to open/read underlying %BAM or
    ///         BAI files.
    ///
    GenomicIntervalQuery(const Data::GenomicInterval& interval, const DataSet& dataset);
    GenomicIntervalQuery(const Data::GenomicInterval& interval, const DataSet& dataset,
                         const BaiIndexCache& cache);

    GenomicIntervalQuery(GenomicIntervalQuery&&) noexcept;
    GenomicIntervalQuery& operator=(GenomicIntervalQuery&&) noexcept;
    ~GenomicIntervalQuery() override;

public:
    /// \brief Main iteration point for record access.
    ///
    /// Most client code should not need to use this method directly. Use
    /// iterators instead.
    ///
    bool GetNext(BamRecord& r) override;

public:
    /// \brief Sets a new genomic interval.
    ///
    /// This allows the same dataset/query to be re-used over multiple regions of
    /// interest:
    ///
    /// \include code/GenomicIntervalQuery_Reuse.txt
    ///
    /// \param[in] interval new genomic interval
    /// \returns reference to this query
    ///
    GenomicIntervalQuery& Interval(const Data::GenomicInterval& interval);

    /// \returns Current genomic interval active on this query.
    const Data::GenomicInterval& Interval() const;

private:
    class GenomicIntervalQueryPrivate;
    std::unique_ptr<GenomicIntervalQueryPrivate> d_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_GENOMICINTERVALQUERY_H
