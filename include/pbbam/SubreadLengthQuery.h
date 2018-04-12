// File Description
/// \file SubreadLengthQuery.h
/// \brief Defines the SubreadLengthQuery class.
//
// Author: Derek Barnett

#ifndef SUBREADLENGTHQUERY_H
#define SUBREADLENGTHQUERY_H

#include <cstdint>
#include <vector>
#include "pbbam/Compare.h"
#include "pbbam/Config.h"
#include "pbbam/internal/QueryBase.h"

namespace PacBio {
namespace BAM {

/// \brief The SubreadLengthQuery class provides iterable access to a DataSet's
///        %BAM records, limiting results to those matching a subread length
///        criterion.
///
/// Example:
/// \include code/SubreadLengthQuery.txt
///
/// \note Currently, all %BAM files must have a corresponding ".pbi" index file.
///       Use BamFile::EnsurePacBioIndexExists before creating the query if one
///       may not be present.
///
class PBBAM_EXPORT SubreadLengthQuery : public internal::IQuery
{
public:
    /// \brief Creates a new SubreadLengthQuery, limiting record results to only
    ///        those matching a subread length criterion.
    ///
    /// \param[in] length       subread length value
    /// \param[in] compareType  compare operator
    /// \param[in] dataset      input data source(s)
    ///
    /// \throws std::runtime_error on failure to open/read underlying %BAM or PBI
    ///         files.
    ///
    SubreadLengthQuery(const int32_t length, const Compare::Type compareType,
                       const DataSet& dataset);

    ~SubreadLengthQuery();

public:
    /// \brief Main iteration point for record access.
    ///
    /// Most client code should not need to use this method directly. Use
    /// iterators instead.
    ///
    bool GetNext(BamRecord& r) override;

    uint32_t NumReads() const;

private:
    struct SubreadLengthQueryPrivate;
    std::unique_ptr<SubreadLengthQueryPrivate> d_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // SUBREADLENGTHQUERY_H
