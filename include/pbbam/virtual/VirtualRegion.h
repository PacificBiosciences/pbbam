// File Description
/// \file VirtualRegion.h
/// \brief Defines the VirtualRegion class.
//
// Author: Armin Töpfer

#ifndef VIRTUALREGION_H
#define VIRTUALREGION_H

#include "pbbam/Config.h"
#include "pbbam/LocalContextFlags.h"
#include "pbbam/virtual/VirtualRegionType.h"

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
    VirtualRegion(const VirtualRegionType type, const int beginPos, const int endPos,
                  const int score = 0);

    /// \brief Creates a virtual region with type/position info, as well as context & barcode.
    ///
    VirtualRegion(const VirtualRegionType type, const int beginPos, const int endPos,
                  const LocalContextFlags cxTag, const int barcodeLeft, const int barcodeRight,
                  const int score = 0);

    VirtualRegion() = default;
    VirtualRegion(const VirtualRegion&) = default;
    VirtualRegion(VirtualRegion&&) = default;
    VirtualRegion& operator=(const VirtualRegion&) = default;  // un-"delete"-ed for SWIG
    VirtualRegion& operator=(VirtualRegion&&) = default;
    ~VirtualRegion() = default;

    bool operator==(const VirtualRegion& v1) const;
};

inline VirtualRegion::VirtualRegion(const VirtualRegionType type, const int beginPos,
                                    const int endPos, const int score)
    : type(type), beginPos(beginPos), endPos(endPos), cxTag(), score(score)
{
}

inline VirtualRegion::VirtualRegion(const VirtualRegionType type, const int beginPos,
                                    const int endPos, const LocalContextFlags cxTag,
                                    const int barcodeLeft, const int barcodeRight, const int score)
    : type(type)
    , beginPos(beginPos)
    , endPos(endPos)
    , cxTag(cxTag)
    , barcodeLeft(barcodeLeft)
    , barcodeRight(barcodeRight)
    , score(score)
{
}

inline bool VirtualRegion::operator==(const VirtualRegion& v1) const
{
    return (v1.type == this->type && v1.beginPos == this->beginPos && v1.endPos == this->endPos);
}

}  // namespace BAM
}  // namespace PacBio

#endif  // VIRTUALREGION_H
