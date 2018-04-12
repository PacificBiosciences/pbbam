// File Description
/// \file ZmwType.h
/// \brief Defines the ZmwType enum.
//
// Author: Armin Töpfer

#ifndef ZMWTYPE_H
#define ZMWTYPE_H

#include "pbbam/Config.h"

namespace PacBio {
namespace BAM {

/// \brief This enum defines the different ZMW categories of scraps
///
enum class ZmwType : char
{
    CONTROL = 'C',
    MALFORMED = 'M',
    NORMAL = 'N',
    SENTINEL = 'S'
};

}  // namespace BAM
}  // namespace PacBio

#endif  // ZMWTYPE_H
