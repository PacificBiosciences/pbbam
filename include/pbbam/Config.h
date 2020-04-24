// File Description
/// \file Config.h
/// \brief Defines library-wide macros & global variables.
//
// Author: Derek Barnett

#ifndef PBBAM_CONFIG_H
#define PBBAM_CONFIG_H

/// Library Import/Export
#ifndef PBBAM_EXPORT
#if defined(WIN32)
#define PBBAM_EXPORT __declspec(dllimport)
#else
#define PBBAM_EXPORT
#endif
#endif

/// Switch for warnings for the pbbam -> pbcopper Data:: move
#ifdef PACBIO_NODEPRECATED_API
#define PBBAM_DEPRECATED [[deprecated("Use the version from pbcopper in Data::")]]
#else
#define PBBAM_DEPRECATED
#endif

/// Switch for warnings for the FrameCodec -> FrameEncodingType move
#ifdef PACBIO_NODEPRECATED_FRAMES
#define PBBAM_DEPRECATED_FRAMES [[deprecated("Use FrameCodec instead.")]]
#else
#define PBBAM_DEPRECATED_FRAMES
#endif

/// Disable use of getrandom(), which requires Linux kernel 3.17+.
/// This define allows use of getentropy() in glibc 2.25+, otherwise
/// fallback to 'posix' provider
#ifndef BOOST_UUID_RANDOM_PROVIDER_DISABLE_GETRANDOM
#define BOOST_UUID_RANDOM_PROVIDER_DISABLE_GETRANDOM
#endif

namespace PacBio {
namespace BAM {

/// Sets the desired verbosity level of htslib warnings.
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

/// \return true if runtime htslib is >= v1.7
bool DoesHtslibSupportLongCigar();

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_CONFIG_H
