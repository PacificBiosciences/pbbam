#include "PbbamInternalConfig.h"

#include "ErrnoReason.h"

#include <cerrno>
#include <cstring>

#include <iostream>

namespace PacBio {
namespace BAM {

void MaybePrintErrnoReason(std::ostream& out)
{
    if (errno != 0) out << "\n  reason: " << std::strerror(errno);
}

}  // namespace BAM
}  // namespace PacBio
