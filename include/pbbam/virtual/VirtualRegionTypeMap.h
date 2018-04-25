// File Description
/// \file VirtualRegionTypeMap.h
/// \brief Defines the VirtualRegionTypeMap class.
//
// Author: Derek Barnett

#ifndef VIRTUALREGIONTYPEMAP_H
#define VIRTUALREGIONTYPEMAP_H

#include <map>

#include "pbbam/Config.h"
#include "pbbam/virtual/VirtualRegionType.h"

namespace PacBio {
namespace BAM {

/// \brief The VirtualRegionTypeMap class provides mapping between char codes and
///        VirtualRegionType enum keys.
///
class VirtualRegionTypeMap
{
public:
    static std::map<char, VirtualRegionType> ParseChar;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // VIRTUALREGIONTYPEMAP_H
