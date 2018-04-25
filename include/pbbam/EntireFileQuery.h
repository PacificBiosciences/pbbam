// File Description
/// \file EntireFileQuery.h
/// \brief Defines the EntireFileQuery class.
//
// Author: Derek Barnett

#ifndef ENTIREFILEQUERY_H
#define ENTIREFILEQUERY_H

#include <memory>
#include "pbbam/internal/QueryBase.h"

namespace PacBio {
namespace BAM {

/// \brief The EntireFileQuery class provides iterable access to a DataSet's
///        %BAM records, reading through the entire contents of each file.
///
/// Input files will be accessed in the order listed in the DataSet.
///
/// \include code/EntireFileQuery.txt
///
/// Iteration is not limited to only 'const' records. The files themselves will
/// not be affected, but individual records may be modified if needed.
///
/// \include code/EntireFileQuery_NonConst.txt
///
/// \note DataSets can be implicitly constructed from %BAM filenames as well.
///       Thus a single %BAM file can be read through using the following:
///
/// \include code/EntireFileQuery_BamFilename.txt
///
class PBBAM_EXPORT EntireFileQuery : public internal::IQuery
{
public:
    /// \brief Creates a new EntireFileQuery, reading through the entire
    ///        contents of a dataset.
    ///
    /// \param[in] dataset  input data source(s)
    /// \throws std::runtime_error on failure to open/read underlying %BAM
    ///         files.
    ///
    EntireFileQuery(const PacBio::BAM::DataSet& dataset);
    ~EntireFileQuery() override;

public:
    /// \brief Main iteration point for record access.
    ///
    /// Most client code should not need to use this method directly. Use
    /// iterators instead.
    ///
    bool GetNext(BamRecord& r) override;

private:
    struct EntireFileQueryPrivate;
    std::unique_ptr<EntireFileQueryPrivate> d_;
};

}  // namespace BAM
}  // namspace PacBio

#endif  // ENTIREFILEQUERY_H
