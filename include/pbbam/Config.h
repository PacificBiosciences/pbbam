// File Description
/// \file Config.h
/// \brief Defines library-wide macros & global variables.
//
// Author: Derek Barnett

#ifndef PBBAM_CONFIG_H
#define PBBAM_CONFIG_H

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

/// \}

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_CONFIG_H
