#ifndef PBBAM_ZMWQUERY_H
#define PBBAM_ZMWQUERY_H

#include <pbbam/Config.h>

#include <pbbam/internal/QueryBase.h>

#include <vector>

#include <cstdint>

namespace PacBio {
namespace BAM {

/// \brief The ZmwQuery class provides iterable access to a DataSet's
///        %BAM records, limiting results to those matching a ZMW hole number
///        whitelist.
///
/// Example:
/// \include code/ZmwQuery.txt
///
/// \note Currently, all %BAM files must have a corresponding ".pbi" index file.
///       Use BamFile::EnsurePacBioIndexExists before creating the query if one
///       may not be present.
///
class PBBAM_EXPORT ZmwQuery : public internal::IQuery
{
public:
    /// \brief Creates a new ZmwQuery, limiting record results to only
    ///        those matching a ZMW hole number criterion.
    ///
    /// \param[in] zmwWhitelist     vector of allowed ZMW hole numbers
    /// \param[in] dataset          input data source(s)
    ///
    /// \throws std::runtime_error on failure to open/read underlying %BAM or
    ///         PBI files.
    ///
    ZmwQuery(std::vector<std::int32_t> zmwWhitelist, const DataSet& dataset);

    ZmwQuery(ZmwQuery&&) noexcept;
    ZmwQuery& operator=(ZmwQuery&&) noexcept;
    ~ZmwQuery() override;

    /// \brief Main iteration point for record access.
    ///
    /// Most client code should not need to use this method directly. Use
    /// iterators instead.
    ///
    bool GetNext(BamRecord& r) override;

private:
    class ZmwQueryPrivate;
    std::unique_ptr<ZmwQueryPrivate> d_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_ZMWQUERY_H
