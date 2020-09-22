#ifndef PBBAM_VIRTUALREGIONTYPEMAP_H
#define PBBAM_VIRTUALREGIONTYPEMAP_H

#include <pbbam/Config.h>

#include <map>

#include <pbbam/virtual/VirtualRegionType.h>

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

#endif  // PBBAM_VIRTUALREGIONTYPEMAP_H
