// Author: Derek Barnett

#ifndef ERRNOREASON_H
#define ERRNOREASON_H

#include "pbbam/Config.h"

#include <iosfwd>

namespace PacBio {
namespace BAM {

void MaybePrintErrnoReason(std::ostream& out);

}  // namespace BAM
}  // namespace PacBio

#endif  // ERRNOREASON_H
