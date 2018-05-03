// File Description
/// \file Accuracy.h
/// \brief Defines the Accuracy class.
//
// Author: Derek Barnett

#ifndef ACCURACY_H
#define ACCURACY_H

#include "pbbam/Config.h"

namespace PacBio {
namespace BAM {

/// \brief The Accuracy class represents the expected accuracy of a BamRecord.
///
/// Values are clamped to fall within [0,1].
///
class PBBAM_EXPORT Accuracy
{
public:
    static const float MIN;  ///< Minimum valid accuracy value [0.0]
    static const float MAX;  ///< Maximum valid accuracy value [1.0]

public:
    /// \name Constructors & Related Methods
    /// \{

    /// Constructs an Accuracy object from a floating-point number.
    ///
    /// \note This is not an \b explicit ctor, to make it as easy as
    ///       possible to use in numeric operations. We really just want
    ///       to make sure that the acceptable range is respected.
    ///
    Accuracy(float accuracy);

    Accuracy(const Accuracy&) = default;
    Accuracy(Accuracy&&) = default;
    Accuracy& operator=(const Accuracy&) = default;
    Accuracy& operator=(Accuracy&&) = default;
    ~Accuracy() = default;

    /// \}

public:
    /// \returns Accuracy as float primitive
    operator float() const;

private:
    float accuracy_;
};

}  // namespace BAM
}  // namespace PacBio

#include "pbbam/internal/Accuracy.inl"

#endif  // ACCURACY_H
