// File Description
/// \file VirtualZmwCompositeReader.h
/// \brief Defines the VirtualZmwCompositeReader class.
//
// Author: Derek Barnett

#ifndef VIRTUALZMWCOMPOSITEREADER_H
#define VIRTUALZMWCOMPOSITEREADER_H

#include <deque>
#include <memory>
#include <string>
#include <utility>

#include "VirtualZmwReader.h"
#include "pbbam/DataSet.h"
#include "pbbam/PbiFilter.h"

namespace PacBio {
namespace BAM {
namespace internal {

/// \brief The VirtualZmwCompositeReader provides an interface for
///        re-stitching "virtual" polymerase reads from their constituent parts,
///        across multiple %BAM resources from a DataSet.
///
/// This class is essentially a DataSet-aware wrapper around
/// VirtualZmwReader, enabling multiple resources as input. See that
/// class's documentation for more info.
///
class PBBAM_EXPORT VirtualZmwCompositeReader
{
public:
    /// \name Constructors & Related Methods
    /// \{

    VirtualZmwCompositeReader(const DataSet& dataset);

    VirtualZmwCompositeReader() = delete;
    VirtualZmwCompositeReader(const VirtualZmwCompositeReader&) = delete;
    VirtualZmwCompositeReader(VirtualZmwCompositeReader&&) = delete;
    VirtualZmwCompositeReader& operator=(const VirtualZmwCompositeReader&) = delete;
    VirtualZmwCompositeReader& operator=(VirtualZmwCompositeReader&&) = delete;
    ~VirtualZmwCompositeReader() = default;

    /// \}

public:
    /// \name Stitched Record Reading
    ///

    /// \returns true if more ZMWs/files are available for reading.
    bool HasNext();

    /// \returns the next stitched polymerase read
    VirtualZmwBamRecord Next();

    /// \returns the next set of reads that belong to one ZMW from one %BAM
    ///          resource (a primary %BAM and/or its scraps file). This enables
    ///          stitching records in a distinct thread.
    ///
    std::vector<BamRecord> NextRaw();

    /// \}

private:
    std::deque<std::pair<std::string, std::string> > sources_;
    std::unique_ptr<VirtualZmwReader> currentReader_;
    PbiFilter filter_;

private:
    void OpenNextReader();
};

}  // namespace internal
}  // namespace BAM
}  // namespace PacBio

#endif  // VIRTUALCOMPOSITEREADER_H
