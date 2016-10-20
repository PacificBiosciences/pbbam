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

// Author: Derek Barnett

#ifndef SEQUENCEUTILS_H
#define SEQUENCEUTILS_H

#include "StringUtils.h"
#include <algorithm>
#include <string>
#include <ctype.h>

namespace PacBio {
namespace BAM {
namespace internal {

inline char Complement(const char character)
{
    static char const complementLookup[] =
    {
        '\0', 'T', 'V', 'G', 'H',
        '\0', '\0', 'C', 'D', '\0',
        '\0', 'M', '\0', 'K', 'N',
        '\0', '\0', '\0', 'Y', 'S',
        'A', 'A', 'B', 'W', '\0', 'R'
    };
    if (character == '-' || character == '*')
        return character;
    return complementLookup[toupper(character) & 0x1f];
}

//inline void Reverse(std::string& s)
//{ std::reverse(s.begin(), s.end()); }

template<typename T>
void Reverse(T& input)
{ std::reverse(input.begin(), input.end()); }

template<typename T>
T MaybeReverse(T&& input, bool reverse)
{ 
    if (reverse) std::reverse(input.begin(), input.end()); 
    return input;
}

template<typename T>
T Reversed(const T& input)
{
    T result = input;
    Reverse(result);
    return result;
}

//inline std::string Reversed(const std::string& input)
//{
//    std::string result = input;
//    Reverse(result);
//    return result;
//}

inline void ReverseComplement(std::string& seq) {

    std::string::iterator sIter = seq.begin();
    std::string::iterator sEnd  = seq.end();
    for ( ; sIter != sEnd; ++sIter )
        *sIter = Complement(*sIter);
    Reverse(seq);
}

inline std::string MaybeReverseComplement(std::string&& seq, bool reverse)
{
    if (reverse) ReverseComplement(seq);
    return seq;
}

/// Reverse complement a DNA sequence case-sensitive
inline void ReverseComplementCaseSens(std::string& seq)
{
    const std::string original = seq;
    int8_t rc_table[128] = {
        4, 4, 4,   4,  4,   4, 4, 4,  4,  4,  4,  4, 4, 4,  4,  4, 4, 4, 4,
        4, 4, 4,   4,  4,   4, 4, 4,  4,  4,  4,  4, 4, 32, 4,  4, 4, 4, 4,
        4, 4, 4,   4,  42,  4, 4, 45, 4,  4,  4,  4, 4, 4,  4,  4, 4, 4, 4,
        4, 4, 4,   4,  4,   4, 4, 4,  84, 4,  71, 4, 4, 4,  67, 4, 4, 4, 4,
        4, 4, 78,  4,  4,   4, 4, 4,  65, 65, 4,  4, 4, 4,  4,  4, 4, 4, 4,
        4, 4, 116, 4,  103, 4, 4, 4,  99, 4,  4,  4, 4, 4,  4,  4, 4, 4, 4,
        4, 4, 97,  97, 4,   4, 4, 4,  4,  4,  4,  4, 4, 4};
    std::string reverseCompl(original.length(), 'N');
    for (uint32_t i = 0; i < original.length(); ++i)
        reverseCompl[original.length()-i-1] = (char)rc_table[(int8_t)original[i]];
    seq = reverseCompl;
}

inline std::string MaybeReverseComplementCaseSens(std::string&& seq, bool reverse)
{
    if (reverse) ReverseComplementCaseSens(seq);
    return seq;
}


inline std::string ReverseComplemented(const std::string& input)
{
    std::string result = input;
    ReverseComplement(result);
    return result;
}

} // namespace internal
} // namespace BAM
} // namespace PacBio

#endif // SEQUENCEUTILS_H
