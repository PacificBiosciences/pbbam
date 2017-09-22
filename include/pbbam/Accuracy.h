// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.
//
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
    static const float MIN; ///< Minimum valid accuracy value [0.0]
    static const float MAX; ///< Maximum valid accuracy value [1.0]

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

    Accuracy(const Accuracy& other) = default;
    Accuracy(Accuracy&& other) = default;
    Accuracy& operator=(const Accuracy& other) = default;
    Accuracy& operator=(Accuracy&& other) = default;
    ~Accuracy() = default;

    /// \}

public:
    /// \returns Accuracy as float primitive
    operator float() const;

private:
    float accuracy_;
};

} // namespace BAM
} // namespace PacBio

#include "pbbam/internal/Accuracy.inl"

#endif // ACCURACY_H
