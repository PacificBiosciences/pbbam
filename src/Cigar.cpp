// Copyright (c) 2014, Pacific Biosciences of California, Inc.
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
#include <sstream>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

// CigarOperation

string CigarOperation::cigarCodes_ = string("MIDNSHP=X");

CigarOperation::CigarOperation(void)
    : type_()
{ }

CigarOperation::CigarOperation(char type, uint32_t length)
    : type_(type)
    , length_(length)
{ }

CigarOperation::CigarOperation(CigarOperationType op, uint32_t length)
    : type_(CigarOperation::TypeToChar(op))
    , length_(length)
{ }

CigarOperation::CigarOperation(const CigarOperation &other)
    : type_(other.type_)
    , length_(other.length_)
{ }

CigarOperation::~CigarOperation(void) { }

char CigarOperation::TypeToChar(const CigarOperationType& type)
{
    return bam_cigar_opchr(static_cast<int>(type));
}

CigarOperationType CigarOperation::CharToType(const char c)
{
    size_t index = cigarCodes_.find(c);
    if (index == string::npos)
        index = 0;
    return static_cast<CigarOperationType>(index);
}

// Cigar

Cigar::Cigar(void)
    : vector<CigarOperation>()
{ }

Cigar::Cigar(const Cigar& other)
    : vector<CigarOperation>(other)
{ }

Cigar::~Cigar(void) { }

Cigar Cigar::FromStdString(const string& stdString)
{
    Cigar result;

    size_t numberStart = 0;
    const size_t numChars = stdString.size();
    for (size_t i = 0; i < numChars; ++i) {
        const char c = stdString.at(i);
        if (!::isdigit(c)) {
            const size_t distance = i - numberStart;
            const uint32_t length = stoul(stdString.substr(numberStart, distance));
            result.push_back(CigarOperation(c, length));
            numberStart = i+1;
        }
    }

    return result;
}

string Cigar::ToStdString(void) const
{
    return Cigar::ToStdString(*this);
}

string Cigar::ToStdString(const Cigar& cigar)
{
    stringstream s;
    auto end  = cigar.cend();
    for (auto iter = cigar.cbegin(); iter != end; ++iter) {
        const CigarOperation& cigar = (*iter);
        s << cigar.Length()
          << cigar.Type();
    }
    return s.str();
}
