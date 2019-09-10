// File Description
/// \file ZmwGroupQuery.h
/// \brief Defines the ZmwGroupQuery class.
//
// Author: Derek Barnett

#ifndef ZMWGROUPQUERY_H
#define ZMWGROUPQUERY_H

#include "pbbam/Config.h"

#include <cstdint>
#include <vector>

#include "pbbam/DataSet.h"
#include "pbbam/internal/QueryBase.h"

namespace PacBio {
namespace BAM {

/// \brief The ZmwGroupQuery class provides iterable access to a DataSet's
///        %BAM records, (optionally limiting results to those matching a ZMW
///        whitelist), grouping those results by hole number.
///
/// Example:
/// \include code/ZmwGroupQuery.txt
///
/// \note %BAM files for whitelisted queries, or with DataSetFilterMode::APPLY,
///       must have a corresponding ".pbi" index file.
///       Use BamFile::EnsurePacBioIndexExists before creating the query if one
///       may not be present.
///
class PBBAM_EXPORT ZmwGroupQuery : public internal::IGroupQuery
{
public:
    ///
    /// \brief Creates a new ZmwGroupQuery that returns BamRecords grouped by ZMW,
    ///        accessing each file in a 'round-robin' fashion.
    ///
    /// \param dataset      input data
    /// \param filterMode   apply/ignore any filters in XML, if present
    ///
    ZmwGroupQuery(const DataSet& dataset,
                  const DataSetFilterMode filterMode = DataSetFilterMode::APPLY);

    /// \brief Creates a new ZmwGroupQuery, limiting record results to only
    ///        those matching a ZMW hole number criterion.
    ///
    /// \param[in] zmwWhitelist     vector of allowed ZMW hole numbers
    /// \param[in] dataset          input data source(s)
    ///
    /// \throws std::runtime_error on failure to open/read underlying %BAM or
    ///         PBI files.
    ///
    ZmwGroupQuery(const std::vector<int32_t>& zmwWhitelist, const DataSet& dataset);
    ~ZmwGroupQuery();

public:
    /// \brief Main iteration point for record access.
    ///
    /// Most client code should not need to use this method directly. Use
    /// iterators instead.
    ///
    bool GetNext(std::vector<BamRecord>& records) override;

private:
    class ZmwGroupQueryPrivate;
    class RoundRobinZmwGroupQuery;
    class WhitelistedZmwGroupQuery;
    std::unique_ptr<ZmwGroupQueryPrivate> d_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // ZMWGROUPQUERY_H
