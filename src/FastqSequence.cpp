// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/FastqSequence.h"

#include <cstdio>
#include <exception>

namespace PacBio {
namespace BAM {

FastqSequence::FastqSequence(std::string name, std::string bases, QualityValues qualities)
    : FastaSequence{std::move(name), std::move(bases)}, qualities_{std::move(qualities)}
{
}

FastqSequence::FastqSequence(std::string name, std::string bases, std::string qualities)
    : FastaSequence{std::move(name), std::move(bases)}
    , qualities_{QualityValues::FromFastq(qualities)}
{
}

const QualityValues& FastqSequence::Qualities() const { return qualities_; }

}  // namespace BAM
}  // namespace PacBio