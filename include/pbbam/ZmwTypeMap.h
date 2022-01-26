#ifndef PBBAM_ZMWTYPEMAP_H
#define PBBAM_ZMWTYPEMAP_H

#include <pbbam/Config.h>

#include <pbbam/ZmwType.h>

#include <map>

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

#endif  // PBBAM_ZMWTYPEMAP_H
