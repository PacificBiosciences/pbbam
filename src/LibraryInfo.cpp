#include "PbbamInternalConfig.h"

#include <pbbam/LibraryInfo.h>

#include "LibraryGitHash.h"
#include "LibraryVersion.h"

#include <pbcopper/LibraryInfo.h>

#include <htslib/hts.h>
#include <zlib.h>

namespace PacBio {
namespace Pbbam {

Library::Info HtslibLibraryInfo() { return {"htslib", std::string{hts_version()}, ""}; }

Library::Info ZlibLibraryInfo() { return {"zlib", std::string{ZLIB_VERSION}, ""}; }

Library::Bundle LibraryBundle()
{
    Library::Bundle bundle{LibraryInfo(), {}};
    bundle += Pbcopper::LibraryBundle();
    bundle += HtslibLibraryInfo();
    bundle += ZlibLibraryInfo();
    return bundle;
}

Library::Info LibraryInfo() { return {"pbbam", Pbbam::ReleaseVersion, Pbbam::LibraryGitSha1}; }

}  // namespace Pbbam
}  // namespace PacBio
