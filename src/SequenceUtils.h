// Author: Derek Barnett

#ifndef SEQUENCEUTILS_H
#define SEQUENCEUTILS_H

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <string>

#include "StringUtils.h"

namespace PacBio {
namespace BAM {
namespace internal {

inline char Complement(const char character)
{
    static char const complementLookup[] = {'\0', 'T',  'V', 'G',  'H', '\0', '\0', 'C',  'D',
                                            '\0', '\0', 'M', '\0', 'K', 'N',  '\0', '\0', '\0',
                                            'Y',  'S',  'A', 'A',  'B', 'W',  '\0', 'R'};
    if (character == '-' || character == '*') return character;
    return complementLookup[toupper(character) & 0x1f];
}

//inline void Reverse(std::string& s)
//{ std::reverse(s.begin(), s.end()); }

template <typename T>
void Reverse(T& input)
{
    std::reverse(input.begin(), input.end());
}

template <typename T>
T MaybeReverse(T&& input, bool reverse)
{
    if (reverse) std::reverse(input.begin(), input.end());
    return input;
}

template <typename T>
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

inline void ReverseComplement(std::string& seq)
{
    std::transform(seq.begin(), seq.end(), seq.begin(), Complement);
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
    constexpr const static int8_t rc_table[128] = {
        4,  4, 4, 4, 4, 4, 4,  4,  4, 4, 4, 4, 4,  4,  4, 4,  4,  4, 4, 4,   4, 4,   4, 4, 4, 4,
        4,  4, 4, 4, 4, 4, 32, 4,  4, 4, 4, 4, 4,  4,  4, 4,  42, 4, 4, 45,  4, 4,   4, 4, 4, 4,
        4,  4, 4, 4, 4, 4, 4,  4,  4, 4, 4, 4, 4,  84, 4, 71, 4,  4, 4, 67,  4, 4,   4, 4, 4, 4,
        78, 4, 4, 4, 4, 4, 65, 65, 4, 4, 4, 4, 4,  4,  4, 4,  4,  4, 4, 116, 4, 103, 4, 4, 4, 99,
        4,  4, 4, 4, 4, 4, 4,  4,  4, 4, 4, 4, 97, 97, 4, 4,  4,  4, 4, 4,   4, 4,   4, 4};
    std::string reverseCompl(original.length(), 'N');
    for (uint32_t i = 0; i < original.length(); ++i)
        reverseCompl[original.length() - i - 1] =
            static_cast<char>(rc_table[static_cast<int8_t>(original[i])]);
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

}  // namespace internal
}  // namespace BAM
}  // namespace PacBio

#endif  // SEQUENCEUTILS_H
