// File Description
/// \file VirtualRegionTypeMap.cpp
/// \brief Implements the VirtualRegionTypeMap class.
//
// Author: Armin Töpfer

#include "PbbamInternalConfig.h"

#include "pbbam/virtual/VirtualRegion.h"

namespace PacBio {
namespace BAM {

VirtualRegion::VirtualRegion(const VirtualRegionType type_, const int beginPos_, const int endPos_,
                             const int score_)
    : type{type_}, beginPos{beginPos_}, endPos{endPos_}, cxTag{}, score{score_}
{
}

VirtualRegion::VirtualRegion(const VirtualRegionType type_, const int beginPos_, const int endPos_,
                             const LocalContextFlags cxTag_, const int barcodeLeft_,
                             const int barcodeRight_, const int score_)
    : type{type_}
    , beginPos{beginPos_}
    , endPos{endPos_}
    , cxTag{cxTag_}
    , barcodeLeft{barcodeLeft_}
    , barcodeRight{barcodeRight_}
    , score{score_}
{
}

VirtualRegion::VirtualRegion() = default;

VirtualRegion::VirtualRegion(const VirtualRegion&) = default;

VirtualRegion::VirtualRegion(VirtualRegion&&) = default;

VirtualRegion& VirtualRegion::operator=(const VirtualRegion&) = default;

VirtualRegion& VirtualRegion::operator=(VirtualRegion&&) = default;

VirtualRegion::~VirtualRegion() = default;

bool VirtualRegion::operator==(const VirtualRegion& v1) const
{
    return (v1.type == this->type && v1.beginPos == this->beginPos && v1.endPos == this->endPos);
}

}  // namespace BAM
}  // namespace PacBio