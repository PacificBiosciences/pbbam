// Author: Derek Barnett

#ifndef PBINDEXDUMP_CPPFORMATTER_H
#define PBINDEXDUMP_CPPFORMATTER_H

#include "PbIndexDumpSettings.h"

namespace PacBio {
namespace PbIndexDump {

struct CppFormatter
{
    static void Run(const Settings& settings);
};

}  // namespace PbIndexDump
}  // namespace PacBio

#endif  // PBINDEXDUMP_CPPFORMATTER_H
