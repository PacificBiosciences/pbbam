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

#include "pbbam/Cigar.h"
#include <htslib/sam.h>
#include <array>
#include <sstream>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

namespace internal {

typedef array<CigarOperationType, 128> CigarLookup;

static
CigarLookup InitCigarLookup(void)
{
    CigarLookup cl;
    cl.fill(CigarOperationType::UNKNOWN);
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

// ----------------
// CigarOperation

CigarOperationType CigarOperation::CharToType(const char c)
{   return (static_cast<uint8_t>(c) >= 128 ? CigarOperationType::UNKNOWN
                                           : internal::cigarLookup_[c] );
}

char CigarOperation::TypeToChar(const CigarOperationType& type)
{ return bam_cigar_opchr(static_cast<int>(type)); }

// ----------------
// Cigar

Cigar::Cigar(const string& cigarString)
    : vector<CigarOperation>()
{
    size_t numberStart = 0;
    const size_t numChars = cigarString.size();
    for (size_t i = 0; i < numChars; ++i) {
        const char c = cigarString.at(i);
        if (!::isdigit(c)) {
            const size_t distance = i - numberStart;
            const uint32_t length = stoul(cigarString.substr(numberStart, distance));
            push_back(CigarOperation(c, length));
            numberStart = i+1;
        }
    }
}

string Cigar::ToStdString(void) const
{
    stringstream s;
    const auto end  = this->cend();
    for (auto iter = this->cbegin(); iter != end; ++iter) {
        const CigarOperation& cigar = (*iter);
        s << cigar.Length()
          << cigar.Char();
    }
    return s.str();
}
