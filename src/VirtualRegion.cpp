// File Description
/// \file VirtualRegionTypeMap.cpp
/// \brief Implements the VirtualRegionTypeMap class.
//
// Author: Armin Töpfer

#include "PbbamInternalConfig.h"

#include <cassert>

#include <tuple>
#include <type_traits>

#include "pbbam/virtual/VirtualRegion.h"

namespace PacBio {
namespace BAM {

static_assert(std::is_copy_constructible<VirtualRegion>::value,
              "VirtualRegion(const VirtualRegion&) is not = default");
static_assert(std::is_copy_assignable<VirtualRegion>::value,
              "VirtualRegion& operator=(const VirtualRegion&) is not = default");

static_assert(std::is_nothrow_move_constructible<VirtualRegion>::value,
              "VirtualRegion(VirtualRegion&&) is not = noexcept");
static_assert(std::is_nothrow_move_assignable<VirtualRegion>::value,
              "VirtualRegion& operator=(VirtualRegion&&) is not = noexcept");

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

bool VirtualRegion::operator==(const VirtualRegion& v1) const
{
    return std::tie(type, beginPos, endPos) == std::tie(v1.type, v1.beginPos, v1.endPos);
}

}  // namespace BAM
}  // namespace PacBio
