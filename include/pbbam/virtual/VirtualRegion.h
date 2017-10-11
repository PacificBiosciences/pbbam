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
/// \file VirtualRegion.h
/// \brief Defines the VirtualRegion class.
//
// Author: Armin Töpfer

#ifndef VIRTUALREGION_H
#define VIRTUALREGION_H

#include "pbbam/Config.h"
#include "pbbam/virtual/VirtualRegionType.h"
#include "pbbam/LocalContextFlags.h"

namespace PacBio {
namespace BAM {

/// \brief The VirtualRegion represents an annotation of a polymerase region.
///
struct VirtualRegion
{
public:
    VirtualRegionType type;
    int beginPos;
    int endPos;
    LocalContextFlags cxTag = LocalContextFlags::NO_LOCAL_CONTEXT;
    int barcodeLeft = -1;
    int barcodeRight = -1;
    int score = 0;

public:
    /// \brief Creates a virtual region with basic type & position info.
    ///
    VirtualRegion(const VirtualRegionType type, 
                  const int beginPos,
                  const int endPos,
                  const int score = 0);

    /// \brief Creates a virtual region with type/position info, as well as context & barcode.
    ///
    VirtualRegion(const VirtualRegionType type, 
                  const int beginPos,
                  const int endPos, 
                  const LocalContextFlags cxTag, 
                  const int barcodeLeft,
                  const int barcodeRight,
                  const int score = 0);

    VirtualRegion() = default;
    VirtualRegion(const VirtualRegion&) = default;
    VirtualRegion(VirtualRegion&&) = default;
    VirtualRegion& operator=(const VirtualRegion&) = default; // un-"delete"-ed for SWIG
    VirtualRegion& operator=(VirtualRegion&&) = default;
    ~VirtualRegion() = default;

    bool operator==(const VirtualRegion &v1) const;

};

inline VirtualRegion::VirtualRegion(const VirtualRegionType type,
                                    const int beginPos,
                                    const int endPos,
                                    const int score)
    : type(type)
    , beginPos(beginPos)
    , endPos(endPos), cxTag()
    , score(score)
{}

inline VirtualRegion::VirtualRegion(const VirtualRegionType type,
                                    const int beginPos,
                                    const int endPos,
                                    const LocalContextFlags cxTag,
                                    const int barcodeLeft,
                                    const int barcodeRight,
                                    const int score)
    : type(type)
    , beginPos(beginPos)
    , endPos(endPos)
    , cxTag(cxTag)
    , barcodeLeft(barcodeLeft)
    , barcodeRight(barcodeRight)
    , score(score)
{}

inline bool VirtualRegion::operator==(const VirtualRegion& v1) const
{
    return (v1.type == this->type &&
            v1.beginPos == this->beginPos &&
            v1.endPos == this->endPos);
}

} // namespace BAM
} // namespace PacBio

#endif // VIRTUALREGION_H
