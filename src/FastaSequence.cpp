// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/FastaSequence.h"

#include <cstdio>

#include <exception>
#include <tuple>

#include <boost/algorithm/string.hpp>

namespace PacBio {
namespace BAM {

FastaSequence::FastaSequence(std::string name, std::string bases)
    : name_{std::move(name)}, bases_{std::move(bases)}
{
    boost::algorithm::trim(name_);
    boost::algorithm::trim(bases_);
}

FastaSequence::FastaSequence() = default;

FastaSequence::FastaSequence(const FastaSequence&) = default;

FastaSequence::FastaSequence(FastaSequence&&) noexcept = default;

FastaSequence& FastaSequence::operator=(const FastaSequence&) = default;

FastaSequence& FastaSequence::operator=(FastaSequence&&) PBBAM_NOEXCEPT_MOVE_ASSIGN = default;

FastaSequence::~FastaSequence() = default;

const std::string& FastaSequence::Bases() const { return bases_; }

const std::string& FastaSequence::Name() const { return name_; }

bool FastaSequence::operator==(const FastaSequence& other) const
{
    return std::tie(name_, bases_) == std::tie(other.name_, other.bases_);
}

bool FastaSequence::operator!=(const FastaSequence& other) const { return !(*this == other); }

}  // namespace BAM
}  // namespace PacBio
