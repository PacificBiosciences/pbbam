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
    VirtualRegion(const VirtualRegionType type_, const int beginPos_, const int endPos_,
                  const int score_ = 0);

    /// \brief Creates a virtual region with type/position info, as well as context & barcode.
    ///
    VirtualRegion(const VirtualRegionType type_, const int beginPos_, const int endPos_,
                  const LocalContextFlags cxTag_, const int barcodeLeft_, const int barcodeRight_,
                  const int score_ = 0);

    VirtualRegion() = default;
    VirtualRegion(const VirtualRegion&) = default;
    VirtualRegion(VirtualRegion&&) = default;
    VirtualRegion& operator=(const VirtualRegion&) = default;  // un-"delete"-ed for SWIG
    VirtualRegion& operator=(VirtualRegion&&) = default;
    ~VirtualRegion() = default;

    bool operator==(const VirtualRegion& v1) const;
};

inline VirtualRegion::VirtualRegion(const VirtualRegionType type_, const int beginPos_,
                                    const int endPos_, const int score_)
    : type{type_}, beginPos{beginPos_}, endPos{endPos_}, cxTag{}, score{score_}
{
}

inline VirtualRegion::VirtualRegion(const VirtualRegionType type_, const int beginPos_,
                                    const int endPos_, const LocalContextFlags cxTag_,
                                    const int barcodeLeft_, const int barcodeRight_,
                                    const int score_)
    : type{type_}
    , beginPos{beginPos_}
    , endPos{endPos_}
    , cxTag{cxTag_}
    , barcodeLeft{barcodeLeft_}
    , barcodeRight{barcodeRight_}
    , score{score_}
{
}

inline bool VirtualRegion::operator==(const VirtualRegion& v1) const
{
    return (v1.type == this->type && v1.beginPos == this->beginPos && v1.endPos == this->endPos);
}

}  // namespace BAM
}  // namespace PacBio

#endif  // VIRTUALREGION_H
