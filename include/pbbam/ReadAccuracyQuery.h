// File Description
/// \file ReadAccuracyQuery.h
/// \brief Defines the ReadAccuracyQuery class.
//
// Author: Derek Barnett

#ifndef READACCURACYQUERY_H
#define READACCURACYQUERY_H

#include <vector>
#include "pbbam/Accuracy.h"
#include "pbbam/Compare.h"
#include "pbbam/Config.h"
#include "pbbam/internal/QueryBase.h"

namespace PacBio {
namespace BAM {

/// \brief The ReadAccuracyQuery class provides iterable access to a DataSet's
///        %BAM records, limiting results to those matching a read accuracy
///        criterion.
///
/// Example:
/// \include code/ReadAccuracyQuery.txt
///
/// \note Currently, all %BAM files must have a corresponding ".pbi" index file.
///       Use BamFile::EnsurePacBioIndexExists before creating the query if one
///       may not be present.
///
class PBBAM_EXPORT ReadAccuracyQuery : public internal::IQuery
{
public:
    /// \brief Creates a new ReadAccuracyQuery, limiting record results to only
    ///        those matching a read accuracy criterion.
    ///
    /// \param[in] accuracy     read accuracy value
    /// \param[in] compareType  compare operator
    /// \param[in] dataset      input data source(s)
    ///
    /// \sa BamRecord::ReadAccuracy
    ///
    /// \throws std::runtime_error on failure to open/read underlying %BAM or PBI
    ///         files.
    ///
    ReadAccuracyQuery(const Accuracy accuracy, const Compare::Type compareType,
                      const DataSet& dataset);

    ~ReadAccuracyQuery() override;

public:
    /// \brief Main iteration point for record access.
    ///
    /// Most client code should not need to use this method directly. Use
    /// iterators instead.
    ///
    bool GetNext(BamRecord& r) override;

    uint32_t NumReads() const;

private:
    struct ReadAccuracyQueryPrivate;
    std::unique_ptr<ReadAccuracyQueryPrivate> d_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // READACCURACYQUERY_H
