// File Description
/// \file MD5.h
/// \brief Defines basic MD5 hash utilities
//
// Author: Brett Bowman

#ifndef MD5_H
#define MD5_H

#include <string>

namespace PacBio {
namespace BAM {

/// \brief MD5 hash of a string as a 32-digit hexadecimal string
///
std::string MD5Hash(const std::string& str);

}  // namespace BAM
}  // namespace PacBio

#endif  // MD5_H
