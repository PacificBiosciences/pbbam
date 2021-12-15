#include "PbbamInternalConfig.h"

#include "ErrnoReason.h"

#include <iostream>

#include <cerrno>
#include <cstring>

namespace PacBio {
namespace BAM {

void MaybePrintErrnoReason(std::ostream& out)
{
    if (errno != 0) {
        out << "\n  reason: " << std::strerror(errno);
    }
}

}  // namespace BAM
}  // namespace PacBio
