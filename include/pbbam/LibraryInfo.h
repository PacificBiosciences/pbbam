// Author: Derek Barnett

#ifndef PBBAM_LIBRARYINFO_H
#define PBBAM_LIBRARYINFO_H

#include "pbbam/Config.h"

#include <pbcopper/library/Bundle.h>
#include <pbcopper/library/Info.h>

namespace PacBio {
namespace Pbbam {

///
/// \return pbbam library info (e.g. name, version)
///
Library::Info LibraryInfo();

///
/// \returns bundle of pbbam library info, plus its dependencies
///
Library::Bundle LibraryBundle();

///
/// \return htslib library info (pbbam dependency)
///
Library::Info HtslibLibraryInfo();

///
/// \return zlib library info (pbbam dependency)
///
Library::Info ZlibLibraryInfo();

}  // namespace Pbbam
}  // namespace PacBio

#endif  // PBBAM_LIBRARYINFO_H
