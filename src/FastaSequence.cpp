#include "PbbamInternalConfig.h"

#include <pbbam/FastaSequence.h>

#include <cstdio>

#include <exception>
#include <tuple>
#include <type_traits>

#include <boost/algorithm/string.hpp>

namespace PacBio {
namespace BAM {

FastaSequence::FastaSequence(std::string name, std::string bases)
    : name_{std::move(name)}, bases_{std::move(bases)}
{
    boost::algorithm::trim(name_);
    boost::algorithm::trim(bases_);
}

const std::string& FastaSequence::Bases() const { return bases_; }

FastaSequence& FastaSequence::Bases(std::string bases)
{
    bases_ = std::move(bases);
    return *this;
}

const std::string& FastaSequence::Name() const { return name_; }

FastaSequence& FastaSequence::Name(std::string name)
{
    name_ = std::move(name);
    return *this;
}

bool FastaSequence::operator==(const FastaSequence& other) const noexcept
{
    return std::tie(name_, bases_) == std::tie(other.name_, other.bases_);
}

bool FastaSequence::operator!=(const FastaSequence& other) const noexcept
{
    return !(*this == other);
}

}  // namespace BAM
}  // namespace PacBio
