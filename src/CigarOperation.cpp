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
/// \file CigarOperation.cpp
/// \brief Implements the CigarOperation class.
//
// Author: Derek Barnett

#include "pbbam/CigarOperation.h"
#include <htslib/sam.h>
#include <array>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

namespace internal {

typedef array<CigarOperationType, 128> CigarLookup;

static
CigarLookup InitCigarLookup(void)
{
    CigarLookup cl;
    cl.fill(CigarOperationType::UNKNOWN_OP);
    cl['M'] = CigarOperationType::ALIGNMENT_MATCH;
    cl['I'] = CigarOperationType::INSERTION;
    cl['D'] = CigarOperationType::DELETION;
    cl['N'] = CigarOperationType::REFERENCE_SKIP;
    cl['S'] = CigarOperationType::SOFT_CLIP;
    cl['H'] = CigarOperationType::HARD_CLIP;
    cl['P'] = CigarOperationType::PADDING;
    cl['='] = CigarOperationType::SEQUENCE_MATCH;
    cl['X'] = CigarOperationType::SEQUENCE_MISMATCH;
    return cl;
}

const static CigarLookup cigarLookup_ = InitCigarLookup();

} // namespace internal

CigarOperationType CigarOperation::CharToType(const char c)
{   return (static_cast<uint8_t>(c) >= 128 ? CigarOperationType::UNKNOWN_OP
                                           : internal::cigarLookup_[c] );
}

char CigarOperation::TypeToChar(const CigarOperationType type)
{ return bam_cigar_opchr(static_cast<int>(type)); }
