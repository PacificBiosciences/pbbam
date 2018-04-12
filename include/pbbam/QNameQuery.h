// File Description
/// \file QNameQuery.h
/// \brief Defines the QNameQuery class.
//
// Author: Derek Barnett

#ifndef QNAMEQUERY_H
#define QNAMEQUERY_H

#include <memory>
#include "pbbam/internal/QueryBase.h"

namespace PacBio {
namespace BAM {

/// \brief The QNameQuery class provides iterable access to a DataSet's records,
///        with each iteration of the query returning a contiguous block of
///        records that share a name.
///
/// There is no random-access here. It is simply a sequential read-through,
/// grouping contiguous results that share a BamRecord::FullName.
///
/// \note The name is not ideal - but for legacy reasons, it will remain as-is
///       for now. It will likely become something more explicit, like
///       "SequentialQNameGroupQuery", so that the name "QNameQuery" will be
///       available for a built-in query on a QNAME filter (or whitelist). This
///       will make it more consistent with other queries (ReadAccuracyQuery,
///       SubreadLengthQuery, ZmwQuery, etc).
///
class PBBAM_EXPORT QNameQuery : public internal::IGroupQuery
{
public:
    /// \brief Creates a new QNameQuery.
    ///
    /// \param[in] dataset      input data source(s)
    ///
    /// \throws std::runtime_error on failure to open/read underlying %BAM files
    ///
    QNameQuery(const DataSet& dataset);
    ~QNameQuery() override;

public:
    /// \brief Main iteration point for record access.
    ///
    /// Most client code should not need to use this method directly. Use
    /// iterators instead.
    ///
    bool GetNext(std::vector<BamRecord>& records) override;

private:
    struct QNameQueryPrivate;
    std::unique_ptr<QNameQueryPrivate> d_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // QNAMEQUERY_H
