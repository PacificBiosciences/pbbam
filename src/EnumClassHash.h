// Copyright (c) 2016, Pacific Biosciences of California, Inc.
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
/// \file EnumClassHash.h
/// \brief Defines the EnumClassHash class.
//
// Author: Derek Barnett

#ifndef ENUMCLASSHASH_H
#define ENUMCLASSHASH_H

#include <cstddef>

namespace PacBio {
namespace BAM {
namespace internal {

///
/// \brief The EnumClassHash struct enables the use of enum class types as keys
///        for std::unordered_map.
///
/// Allows something like:
///
/// \code{.cpp}
///    std::unordered_map<Key_t, Value_t, EnumClassHash> myLookup;
/// \endcode
///
/// where Key_t is an enum class. Without this sort of extra hand-holding to
/// provide a 'manual' hash value, enum classes as keys will fail to compile.
///
/// \note This approach might be unnecessary in C++14, if I understand some of
/// the changes correctly. But this works for C++11 and should continue beyond.
///
/// \sa http://stackoverflow.com/questions/18837857/cant-use-enum-class-as-unordered-map-key
///
struct EnumClassHash
{
    // *** NOTE ***
    //
    // Remove this when we integrate pbcopper.
    // This is a duplicate of pbcopper/utility/EnumClassHash.h
    //

    template<typename T> size_t operator()(const T t) const
    { return static_cast<size_t>(t); }
};

} // namespace internal
} // namespace BAM
} // namespace PacBio

#endif // ENUMCLASSHASH_H
