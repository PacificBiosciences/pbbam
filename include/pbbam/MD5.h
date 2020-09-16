#ifndef PBBAM_MD5_H
#define PBBAM_MD5_H

#include <pbbam/Config.h>

#include <string>

namespace PacBio {
namespace BAM {

/// \brief MD5 hash of a string as a 32-digit hexadecimal string
///
std::string MD5Hash(const std::string& str);

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_MD5_H
