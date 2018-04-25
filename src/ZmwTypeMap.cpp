// File Description
/// \file ZmwTypeMap.cpp
/// \brief Implements the ZmwTypeMap class.
//
// Author: Armin Töpfer

#include "PbbamInternalConfig.h"

#include "pbbam/ZmwTypeMap.h"

namespace PacBio {
namespace BAM {

// clang-format off
std::map<char, ZmwType> ZmwTypeMap::ParseChar
{
    { 'C' , ZmwType::CONTROL   },
    { 'M' , ZmwType::MALFORMED },
    { 'N' , ZmwType::NORMAL    },
    { 'S' , ZmwType::SENTINEL  }
};
// clang-format on

}  // namespace BAM
}  // namespace PacBio
