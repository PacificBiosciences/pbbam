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

// Author: Armin Töpfer

#ifndef POLYMERASEBAMRECORD_H
#define POLYMERASEBAMRECORD_H

#include <vector>
#include <sstream>

#include "pbbam/BamHeader.h"
#include "pbbam/BamRecord.h"
#include "pbbam/Config.h"
#include "pbbam/virtual/VirtualRegion.h"
#include "pbbam/virtual/VirtualRegionType.h"

namespace PacBio {
namespace BAM {

/// This class represents a polymerase read stitched on the fly
/// from subreads|hqregion+scraps.
class VirtualPolymeraseBamRecord : public BamRecord
{
public:
    VirtualPolymeraseBamRecord(std::vector<BamRecord>&& unorderedSources,
                               const BamHeader& header);

    VirtualPolymeraseBamRecord() = delete;
    // Move constructor
    VirtualPolymeraseBamRecord(VirtualPolymeraseBamRecord&&) = default;
    // Copy constructor
    VirtualPolymeraseBamRecord(const VirtualPolymeraseBamRecord&) = default; // un-"delete"-ed for SWIG
    // Move assignment operator
    VirtualPolymeraseBamRecord& operator=(VirtualPolymeraseBamRecord&&) = default;
    // Copy assignment operator
    VirtualPolymeraseBamRecord& operator=(const VirtualPolymeraseBamRecord&) = delete;
    // Destructor
    virtual ~VirtualPolymeraseBamRecord() = default;

public:
    /// Provides bool if a given VirtualRegionType has been annotated
    bool HasVirtualRegionType(const VirtualRegionType regionType) const
    { return virtualRegionsMap_.find(regionType) != virtualRegionsMap_.end(); }

    /// Provides annotations of the polymerase read for a given VirtualRegionType
    std::vector<VirtualRegion> VirtualRegionsTable(const VirtualRegionType regionType) const
    { return virtualRegionsMap_.at(regionType); }

    /// Provides all annotations of the polymerase read as a map
    std::map<VirtualRegionType, std::vector<VirtualRegion>> VirtualRegionsMap() const
    { return virtualRegionsMap_; }

public: // New BamRecord functionality.
    Frames IPDV1Frames(Orientation orientation = Orientation::NATIVE) const;

private:
    std::vector<BamRecord> sources_;
    std::map<VirtualRegionType, std::vector<VirtualRegion>> virtualRegionsMap_;

private:
    void StitchSources();

    /// \brief Appends content of src vector to dst vector using move semantics.
    /// \param[in] src Input vector that will be empty after execution
    /// \param[in,out] dest Output vector that will be appended to
    template <typename T>
    inline void MoveAppend(std::vector<T>& src, std::vector<T>& dst) noexcept
    {
        if (dst.empty())
        {
            dst = std::move(src);
        }
        else
        {
            dst.reserve(dst.size() + src.size());
            std::move(src.begin(), src.end(), std::back_inserter(dst));
            src.clear();
        }
    }

    /// \brief Appends content of src vector to dst vector using move semantics.
    /// \param[in] src Input vector via perfect forwarding
    /// \param[in,out] dest Output vector that will be appended to
    template <typename T>
    inline void MoveAppend(std::vector<T>&& src, std::vector<T>& dst) noexcept
    {
        if (dst.empty())
        {
            dst = std::move(src);
        }
        else
        {
            dst.reserve(dst.size() + src.size());
            std::move(src.begin(), src.end(), std::back_inserter(dst));
            src.clear();
        }
    }
};

} // namespace BAM
} // namespace PacBio

#endif // POLYMERASEBAMRECORD_H
