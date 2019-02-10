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
class VirtualRegion
{
public:
    VirtualRegionType type;
    int beginPos;
    int endPos;
    LocalContextFlags cxTag = LocalContextFlags::NO_LOCAL_CONTEXT;
    int barcodeLeft = -1;
    int barcodeRight = -1;
    int score = 0;

    /// \brief Creates a virtual region with basic type & position info.
    ///
    VirtualRegion(const VirtualRegionType type_, const int beginPos_, const int endPos_,
                  const int score_ = 0);

    /// \brief Creates a virtual region with type/position info, as well as context & barcode.
    ///
    VirtualRegion(const VirtualRegionType type_, const int beginPos_, const int endPos_,
                  const LocalContextFlags cxTag_, const int barcodeLeft_, const int barcodeRight_,
                  const int score_ = 0);

    VirtualRegion();
    VirtualRegion(const VirtualRegion&);
    VirtualRegion(VirtualRegion&&);
    VirtualRegion& operator=(const VirtualRegion&);
    VirtualRegion& operator=(VirtualRegion&&);
    ~VirtualRegion();

    bool operator==(const VirtualRegion& v1) const;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // VIRTUALREGION_H
