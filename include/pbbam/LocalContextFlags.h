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
/// \file LocalContextFlags.h
/// \brief Defines the LocalContextFlags enum & helper method(s).
//
// Author: Lance Hepler

#ifndef LOCALCONTEXTFLAGS_H
#define LOCALCONTEXTFLAGS_H

#include <cstdint>

#include "pbbam/Config.h"

namespace PacBio {
namespace BAM {

/// \brief The LocalContextFlags enum defines the flags that can be used
///        to describe a subread's "local context", i.e. whether it is
///        flanked by barcodes/adapters or its pass orientation.
///
enum LocalContextFlags : uint8_t
{
    NO_LOCAL_CONTEXT = 0,   ///< No context information available
    ADAPTER_BEFORE   = 1,   ///< Adapter precedes subread
    ADAPTER_AFTER    = 2,   ///< Adapter follows subread
    BARCODE_BEFORE   = 4,   ///< Barcode precedes subread
    BARCODE_AFTER    = 8,   ///< Barcode follows subread
    FORWARD_PASS     = 16,  ///< Subread's orientation is 'forward pass'
    REVERSE_PASS     = 32   ///< Subread's orientation is 'reverse pass'
};


/// \returns a LocalContextFlags value containing the result of the bitwise-OR
///          operation of \p lhs and \p rhs.
// constexpr is implicitly inline
constexpr LocalContextFlags operator|(const LocalContextFlags lhs, const LocalContextFlags rhs)
{
    return static_cast<LocalContextFlags>(static_cast<int>(lhs) | static_cast<int>(rhs));
}

} // namespace BAM
} // namespace PacBio

#endif // LOCALCONTEXTFLAGS_H
