// Author: Derek Barnett

#ifndef PBBAM_CCS_CCSHEADER_H
#define PBBAM_CCS_CCSHEADER_H

#include "pbbam/Config.h"

#include <string>

namespace PacBio {
namespace CCS {

struct CCSHeader
{
    std::string MovieName;
    std::string BindingKit;
    std::string SequencingKit;
    std::string BasecallerVersion;
    std::string FrameRate;
};

}  // namespace CCS
}  // namespace PacBio

#endif  // PBBAM_CCS_CCSHEADER_H
