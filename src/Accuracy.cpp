// File Description
/// \file Accuracy.cpp
/// \brief Implements the Accuracy class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/Accuracy.h"

namespace PacBio {
namespace BAM {

const float Accuracy::MIN = 0.0f;
const float Accuracy::MAX = 1.0f;

Accuracy::Accuracy(float accuracy)
{
    if (accuracy < Accuracy::MIN)
        accuracy = Accuracy::MIN;
    else if (accuracy > Accuracy::MAX)
        accuracy = Accuracy::MAX;
    accuracy_ = accuracy;
}

Accuracy::Accuracy(const Accuracy&) = default;

Accuracy::Accuracy(Accuracy&&) noexcept = default;

Accuracy& Accuracy::operator=(const Accuracy&) = default;

Accuracy& Accuracy::operator=(Accuracy&&) noexcept = default;

Accuracy::~Accuracy() = default;

Accuracy::operator float() const { return accuracy_; }

}  // namespace BAM
}  // namespace PacBio
