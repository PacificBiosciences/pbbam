// Copyright (c) 2015, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.
//
// File Description
/// \file VirtualZmwBamRecord.h
/// \brief Defines the VirtualZmwBamRecord class.
//
// Author: Armin Töpfer

#ifndef VirtualZmwBAMRECORD_H
#define VirtualZmwBAMRECORD_H

#include <vector>
#include <sstream>

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
    VirtualZmwBamRecord(std::vector<BamRecord>&& unorderedSources,
                        const BamHeader& header);

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

} // namespace BAM
} // namespace PacBio

#endif // VirtualZmwBAMRECORD_H
