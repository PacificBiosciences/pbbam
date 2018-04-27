// Author: Derek Barnett

#ifndef BAM2SAM_H
#define BAM2SAM_H

#include "Settings.h"

namespace bam2sam {

class PbBam2Sam
{
public:
    static void Run(const Settings& settings);
};

}  // namespace bam2sam

#endif  // PBIBAM2SAM_H
