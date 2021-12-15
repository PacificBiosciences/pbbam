#ifndef PBBAM_VIRTUALREGION_H
#define PBBAM_VIRTUALREGION_H

#include <pbbam/Config.h>

#include <pbbam/LocalContextFlags.h>
#include <pbbam/virtual/VirtualRegionType.h>

namespace PacBio {
namespace BAM {

/// \brief The VirtualRegion represents an annotation of a polymerase region.
///
class VirtualRegion
{
public:
    VirtualRegionType type;
    int beginPos;
    int endPos;
    Data::LocalContextFlags cxTag = Data::LocalContextFlags::NO_LOCAL_CONTEXT;
    int barcodeLeft = -1;
    int barcodeRight = -1;
    int score = 0;

    /// \brief Creates a virtual region with basic type & position info.
    ///
    VirtualRegion(VirtualRegionType type, int beginPos, int endPos, int score = 0);

    /// \brief Creates a virtual region with type/position info, as well as context & barcode.
    ///
    VirtualRegion(VirtualRegionType type, int beginPos, int endPos, Data::LocalContextFlags cxTag,
                  int barcodeLeft, int barcodeRight, int score = 0);

    VirtualRegion() = default;

    bool operator==(const VirtualRegion& v1) const noexcept;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_VIRTUALREGION_H
