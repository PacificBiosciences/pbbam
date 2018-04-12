// File Description
/// \file Accuracy.inl
/// \brief Inline implementations for the Accuracy class.
//
// Author: Derek Barnett

#include "pbbam/Accuracy.h"

namespace PacBio {
namespace BAM {

inline Accuracy::Accuracy(float accuracy)
{
    if (accuracy < Accuracy::MIN)
        accuracy = Accuracy::MIN;
    else if (accuracy > Accuracy::MAX)
        accuracy = Accuracy::MAX;
    accuracy_ = accuracy;
}

inline Accuracy::operator float() const
{ return accuracy_; }

} // namespace BAM
} // namespace PacBio
