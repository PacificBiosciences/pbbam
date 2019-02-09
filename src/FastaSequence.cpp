// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/FastaSequence.h"

#include <cstdio>
#include <exception>

namespace PacBio {
namespace BAM {

FastaSequence::FastaSequence(std::string name, std::string bases)
    : name_{std::move(name)}, bases_{std::move(bases)}
{
    boost::algorithm::trim(name_);
    boost::algorithm::trim(bases_);
}

const std::string& FastaSequence::Bases() const { return bases_; }

const std::string& FastaSequence::Name() const { return name_; }

}  // namespace BAM
}  // namespace PacBio