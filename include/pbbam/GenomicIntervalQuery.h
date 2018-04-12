// File Description
/// \file GenomicIntervalQuery.h
/// \brief Defines the GenomicIntervalQuery class.
//
// Author: Derek Barnett

#ifndef GENOMICINTERVALQUERY_H
#define GENOMICINTERVALQUERY_H

#include <memory>
#include "pbbam/GenomicInterval.h"
#include "pbbam/internal/QueryBase.h"

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
    /// \brief Constructs a new GenomiIntervalQuery, limiting record results to
    ///        only those overalpping a GenomicInterval.
    ///
    /// \param[in] interval genomic interval of interest
    /// \param[in] dataset  input data source(s)
    ///
    /// \throws std::runtime_error on failure to open/read underlying %BAM or
    ///         BAI files.
    ///
    GenomicIntervalQuery(const GenomicInterval& interval, const PacBio::BAM::DataSet& dataset);
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
    GenomicIntervalQuery& Interval(const GenomicInterval& interval);

    /// \returns Current genomic interval active on this query.
    const GenomicInterval& Interval() const;

private:
    struct GenomicIntervalQueryPrivate;
    std::unique_ptr<GenomicIntervalQueryPrivate> d_;
};

}  // namespace BAM
}  // namspace PacBio

#endif  // GENOMICINTERVALQUERY_H
