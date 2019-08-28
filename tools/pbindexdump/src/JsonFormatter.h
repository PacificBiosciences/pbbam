// Author: Derek Barnett

#ifndef PBINDEXDUMP_JSONFORMATTER_H
#define PBINDEXDUMP_JSONFORMATTER_H

#include "PbIndexDumpSettings.h"

namespace PacBio {
namespace PbIndexDump {

struct JsonFormatter
{
    static void Run(const Settings& settings);
};

}  // namespace PbIndexDump
}  // namespace PacBio

#endif  // PBINDEXDUMP_JSONFORMATTER_H
