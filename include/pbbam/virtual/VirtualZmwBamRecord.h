// File Description
/// \file VirtualZmwBamRecord.h
/// \brief Defines the VirtualZmwBamRecord class.
//
// Author: Armin Töpfer

#ifndef VirtualZmwBAMRECORD_H
#define VirtualZmwBAMRECORD_H

#include <sstream>
#include <vector>

#include "pbbam/BamHeader.h"
#include "pbbam/BamRecord.h"
#include "pbbam/Config.h"
#include "pbbam/virtual/VirtualRegion.h"
#include "pbbam/virtual/VirtualRegionType.h"

namespace PacBio {
namespace BAM {

/// \brief The VirtualZmwBamRecord class represents a ZMW read stitched
///        on-the-fly from subreads|hqregion + scraps.
///
class VirtualZmwBamRecord : public BamRecord
{
public:
    /// \name Constructors & Related Methods
    /// \{

    /// \brief Creates a "virtual" ZMW %BAM record, by re-stitching its
    ///        constituent segments.
    ///
    /// \param[in] unorderedSources source data (subreads, scraps, etc.)
    /// \param[in] header           %BAM header to associate with the new record
    ///
    /// \throws std::runtime_error on failure to stitch virtual record
    ///
    VirtualZmwBamRecord(std::vector<BamRecord> unorderedSources, const BamHeader& header);

    VirtualZmwBamRecord() = delete;
    VirtualZmwBamRecord(const VirtualZmwBamRecord&) = default;
    VirtualZmwBamRecord(VirtualZmwBamRecord&&) = default;
    VirtualZmwBamRecord& operator=(const VirtualZmwBamRecord&) = default;
    VirtualZmwBamRecord& operator=(VirtualZmwBamRecord&&) = default;
    virtual ~VirtualZmwBamRecord() = default;

    /// \}

public:
    /// \name Virtual Record Attributes
    ///

    /// \returns true if requested VirtualRegionType has been annotated.
    ///
    bool HasVirtualRegionType(const VirtualRegionType regionType) const;

    /// \returns IPD frame data
    ///
    Frames IPDV1Frames(Orientation orientation = Orientation::NATIVE) const;

    /// \brief Provides all annotations of the polymerase read as a map (type => regions)
    ///
    std::map<VirtualRegionType, std::vector<VirtualRegion>> VirtualRegionsMap() const;

    /// \brief Provides annotations of the polymerase read for a given VirtualRegionType.
    ///
    /// \param[in] regionType  requested region type
    /// \returns regions that match the requested type (empty vector if none found).
    ///
    std::vector<VirtualRegion> VirtualRegionsTable(const VirtualRegionType regionType) const;

    /// \}

private:
    std::vector<BamRecord> sources_;
    std::map<VirtualRegionType, std::vector<VirtualRegion>> virtualRegionsMap_;

private:
    void StitchSources();
};

}  // namespace BAM
}  // namespace PacBio

#endif  // VirtualZmwBAMRECORD_H
