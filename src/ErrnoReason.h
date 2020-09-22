#ifndef PBBAM_ERRNOREASON_H
#define PBBAM_ERRNOREASON_H

#include <pbbam/Config.h>

#include <iosfwd>

namespace PacBio {
namespace BAM {

void MaybePrintErrnoReason(std::ostream& out);

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_ERRNOREASON_H
