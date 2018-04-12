// File Description
/// \file BarcodeQuery.h
/// \brief Defines the BarcodeQuery class.
//
// Author: Derek Barnett

#ifndef BARCODEQUERY_H
#define BARCODEQUERY_H

#include <cstdint>
#include <vector>
#include "pbbam/Config.h"
#include "pbbam/internal/QueryBase.h"

namespace PacBio {
namespace BAM {

/// \brief The BarcodeQuery class provides iterable access to a DataSet's %BAM
///        records, limiting results to those matching a particular barcode.
///
/// Example:
/// \include code/BarcodeQuery.txt
///
/// \note Currently, all %BAM files must have a corresponding ".pbi" index file.
///       Use BamFile::EnsurePacBioIndexExists before creating the query if one
///       may not be present.
///
class PBBAM_EXPORT BarcodeQuery : public internal::IQuery
{
public:
    /// \brief Creates a new BarcodeQuery, limiting record results to only those
    ///        annotated with a particular barcode ID.
    ///
    /// \param[in] barcode  filtering criteria
    /// \param[in] dataset  input data source(s)
    ///
    /// \sa BamRecord::Barcodes
    ///
    /// \throws std::runtime_error on failure to open/read underlying %BAM or PBI
    ///         files.
    ///
    BarcodeQuery(const int16_t barcode, const DataSet& dataset);

    ~BarcodeQuery() override;

public:
    /// \brief Main iteration point for record access.
    ///
    /// Most client code should not need to use this method directly. Use
    /// iterators instead.
    ///
    bool GetNext(BamRecord& r) override;

private:
    struct BarcodeQueryPrivate;
    std::unique_ptr<BarcodeQueryPrivate> d_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // BARCODEQUERY_H
