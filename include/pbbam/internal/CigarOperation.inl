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
/// \file CigarOperation.inl
/// \brief Inline implemenations for the CigarOperation class.
//
// Author: Derek Barnett

#include "pbbam/CigarOperation.h"

namespace PacBio {
namespace BAM {

inline CigarOperation::CigarOperation(char c, uint32_t length)
    : type_(CigarOperation::CharToType(c))
    , length_(length)
{
    #ifndef PBBAM_PERMISSIVE_CIGAR
        if (type_ == CigarOperationType::ALIGNMENT_MATCH)
            throw std::runtime_error("CIGAR operation 'M' is not allowed in PacBio BAM files. Use 'X/=' instead.");
    #endif
}

inline CigarOperation::CigarOperation(CigarOperationType op, uint32_t length)
    : type_(op)
    , length_(length)
{
    #ifndef PBBAM_PERMISSIVE_CIGAR
        if (type_ == CigarOperationType::ALIGNMENT_MATCH)
            throw std::runtime_error("CIGAR operation 'M' is not allowed in PacBio BAM files. Use 'X/=' instead.");
    #endif
}

inline uint32_t CigarOperation::Length() const
{ return length_; }

inline CigarOperation& CigarOperation::Length(const uint32_t length)
{ length_ = length; return *this; }

inline CigarOperationType CigarOperation::Type() const
{ return type_; }

inline CigarOperation &CigarOperation::Type(const CigarOperationType opType)
{ type_ = opType; return *this; }

inline char CigarOperation::Char() const
{ return CigarOperation::TypeToChar(type_); }

inline CigarOperation &CigarOperation::Char(const char opChar)
{ type_ = CigarOperation::CharToType(opChar);return *this; }

inline bool CigarOperation::operator==(const CigarOperation& other) const
{ return type_ == other.type_ && length_ == other.length_; }

inline bool CigarOperation::operator!=(const CigarOperation& other) const
{ return !(*this == other); }

} // namespace BAM
} // namespace PacBio
