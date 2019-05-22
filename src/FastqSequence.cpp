// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/FastqSequence.h"

#include <cstdio>
#include <exception>
#include <tuple>

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

FastqSequence::FastqSequence() = default;

FastqSequence::FastqSequence(const FastqSequence&) = default;

FastqSequence::FastqSequence(FastqSequence&&) noexcept = default;

FastqSequence& FastqSequence::operator=(const FastqSequence&) = default;

FastqSequence& FastqSequence::operator=(FastqSequence&&) noexcept(
    std::is_nothrow_move_assignable<FastaSequence>::value) = default;

FastqSequence::~FastqSequence() = default;

const QualityValues& FastqSequence::Qualities() const { return qualities_; }

FastqSequence& FastqSequence::Qualities(QualityValues quals)
{
    qualities_ = std::move(quals);
    return *this;
}

bool FastqSequence::operator==(const FastqSequence& other) const
{
    return std::tie(Name(), Bases(), qualities_) ==
           std::tie(other.Name(), other.Bases(), other.qualities_);
}

bool FastqSequence::operator!=(const FastqSequence& other) const { return !(*this == other); }

}  // namespace BAM
}  // namespace PacBio
