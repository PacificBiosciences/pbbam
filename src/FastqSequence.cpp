// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/FastqSequence.h"

#include <cassert>
#include <cstdio>

#include <numeric>
#include <stdexcept>
#include <tuple>
#include <type_traits>

namespace PacBio {
namespace BAM {

static_assert(std::is_copy_constructible<FastqSequence>::value,
              "FastqSequence(const FastqSequence&) is not = default");
static_assert(std::is_copy_assignable<FastqSequence>::value,
              "FastqSequence& operator=(const FastqSequence&) is not = default");

static_assert(std::is_nothrow_move_constructible<FastqSequence>::value,
              "FastqSequence(FastqSequence&&) is not = noexcept");
static_assert(std::is_nothrow_move_assignable<FastqSequence>::value ==
                  std::is_nothrow_move_assignable<FastaSequence>::value,
              "");

FastqSequence::FastqSequence(std::string name, std::string bases, Data::QualityValues qualities)
    : FastaSequence{std::move(name), std::move(bases)}, qualities_{std::move(qualities)}
{
}

FastqSequence::FastqSequence(std::string name, std::string bases, std::string qualities)
    : FastaSequence{std::move(name), std::move(bases)}
    , qualities_{Data::QualityValues::FromFastq(qualities)}
{
}

float FastqSequence::AverageBaseQuality() const
{
    if (qualities_.empty()) {
        throw std::runtime_error{
            "[pbbam] FASTQ sequence ERROR: cannot calculate average base quality from "
            "empty sequence"};
    }

    const float total = std::accumulate(
        qualities_.cbegin(), qualities_.cend(), 0.0f,
        [](float currentTotal, const Data::QualityValue& qv) { return currentTotal + qv; });
    return total / qualities_.size();
}

const Data::QualityValues& FastqSequence::Qualities() const { return qualities_; }

FastqSequence& FastqSequence::Qualities(QualityValues quals)
{
    qualities_ = std::move(quals);
    return *this;
}

bool FastqSequence::operator==(const FastqSequence& other) const noexcept
{
    return std::tie(Name(), Bases(), qualities_) ==
           std::tie(other.Name(), other.Bases(), other.qualities_);
}

bool FastqSequence::operator!=(const FastqSequence& other) const noexcept
{
    return !(*this == other);
}

}  // namespace BAM
}  // namespace PacBio
