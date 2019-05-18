// File Description
/// \file Config.h
/// \brief Defines library-wide macros & global variables.
//
// Author: Derek Barnett

#ifndef PBBAM_CONFIG_H
#define PBBAM_CONFIG_H

#include <string>
#include <type_traits>

// string& operator= (string&& str)
// has been made noexcept only in LWG 2063:
//   http://cplusplus.github.io/LWG/lwg-defects.html#2063
// GCC 6+'s "C++11"/SSO std::string is noexcept, but the COW std::string used
// by RHEL devtoolset-6 isn't, and as such we have to make this conditional
// See also
//   https://gcc.gnu.org/bugzilla/show_bug.cgi?id=58265#c13
#define PBBAM_NOEXCEPT_MOVE_ASSIGN noexcept(std::is_nothrow_move_assignable<std::string>::value)

/// \name Library Import/Export
/// \{

#ifndef PBBAM_EXPORT
#if defined(WIN32)
#define PBBAM_EXPORT __declspec(dllimport)
#else
#define PBBAM_EXPORT
#endif
#endif

/// \}

namespace PacBio {
namespace BAM {

/// \name Verbosity Settings
/// \{

/// \brief Sets the desired verbosity level of htslib warnings.
///
/// Change this value to allow debug/warning statements from htslib itself.
/// The valid range seems to be [0-3], where 0 indicates OFF, and 3 is the
/// most verbose.
///
/// By default, pbbam disables htslib statements to keep output channels clean.
/// We rely on exceptions & their associated messages instead.
///
/// This global variable is obviously not thread-safe by any means. But as a
/// debug flag, it is unlikely to cause any real issues. The worst case would be
/// unexpected presence/absence of output statements.
///
extern int HtslibVerbosity;

///
/// \brief DoesHtslibSupportLongCigar
///
/// \return true if runtime htslib is >= v1.7
///
bool DoesHtslibSupportLongCigar();

/// \}

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_CONFIG_H
