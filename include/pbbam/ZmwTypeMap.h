// File Description
/// \file ZmwTypeMap.h
/// \brief Defines the ZmwTypeMap class.
//
// Author: Armin Töpfer

#ifndef ZMWTYPEMAP_H
#define ZMWTYPEMAP_H

#include <map>

#include "pbbam/Config.h"
#include "pbbam/ZmwType.h"

namespace PacBio {
namespace BAM {

/// \brief The ZmwTypeMap class provides mapping between char codes and
///        ZmwType enum keys.
///
class ZmwTypeMap
{
public:
    static std::map<char, ZmwType> ParseChar;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // ZMWTYPEMAP_H
