#include "PbbamInternalConfig.h"

#include <pbbam/virtual/VirtualRegion.h>

#include <tuple>
#include <type_traits>

namespace PacBio {
namespace BAM {

VirtualRegion::VirtualRegion(const VirtualRegionType type_, const int beginPos_, const int endPos_,
                             const int score_)
    : type{type_}, beginPos{beginPos_}, endPos{endPos_}, cxTag{}, score{score_}
{}

VirtualRegion::VirtualRegion(const VirtualRegionType type_, const int beginPos_, const int endPos_,
                             const Data::LocalContextFlags cxTag_, const int barcodeLeft_,
                             const int barcodeRight_, const int score_)
    : type{type_}
    , beginPos{beginPos_}
    , endPos{endPos_}
    , cxTag{cxTag_}
    , barcodeLeft{barcodeLeft_}
    , barcodeRight{barcodeRight_}
    , score{score_}
{}

bool VirtualRegion::operator==(const VirtualRegion& v1) const noexcept
{
    return std::tie(type, beginPos, endPos) == std::tie(v1.type, v1.beginPos, v1.endPos);
}

}  // namespace BAM
}  // namespace PacBio
