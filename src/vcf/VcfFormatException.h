#ifndef PBBAM_VCF_VCFFORMATEXCEPTION_H
#define PBBAM_VCF_VCFFORMATEXCEPTION_H

#include "pbbam/Config.h"

#include <stdexcept>

namespace PacBio {
namespace VCF {

class VcfFormatException : public std::runtime_error
{
public:
    VcfFormatException(std::string reason)
        : std::runtime_error{"[pbbam] VCF format ERROR: " + std::move(reason)}
    {}
};

}  // namespace VCF
}  // namespace PacBio

#endif  // PBBAM_VCF_VCFFORMATEXCEPTION_H
