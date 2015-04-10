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

#ifndef CIGAR_H
#define CIGAR_H

#include "pbbam/CigarOperation.h"
#include "pbbam/Config.h"
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

class PBBAM_EXPORT Cigar : public std::vector<CigarOperation>
{

public:
    /// \name Static Constructor
    /// \{

    /// Creates a Cigar object from SAM/BAM string input
    ///
    /// \param [in] stdString SAM/BAM formatted CIGAR data
    /// \returns Cigar object representing the input data
    static Cigar FromStdString(const std::string& stdString);

    /// \}

public:
    /// \name Constructors & Related Methods
    /// \{

    Cigar(void);
    Cigar(const std::string& cigarString);
    Cigar(const Cigar& other);
    Cigar(Cigar&& other);
    Cigar& operator=(const Cigar& other);
    Cigar& operator=(Cigar&& other);
    ~Cigar(void);

    /// \}

public:
    /// \name Conversion Methods
    /// \{

    /// Converts Cigar object data to SAM/BAM formatted string
    ///
    /// \returns SAM/BAM formatted std::string
    std::string ToStdString(void) const;

    /// \}
};

inline Cigar::Cigar(void)
    : std::vector<CigarOperation>()
{ }

inline Cigar::Cigar(const Cigar& other)
    : std::vector<CigarOperation>(other)
{ }

inline Cigar::Cigar(Cigar&& other)
    : std::vector<CigarOperation>(std::move(other))
{ }

inline Cigar& Cigar::operator=(const Cigar& other)
{ *this = other; return *this; }

inline Cigar& Cigar::operator=(Cigar&& other)
{ *this = std::move(other); return *this; }

inline Cigar::~Cigar(void) { }

inline Cigar Cigar::FromStdString(const std::string& stdString)
{ return Cigar(stdString); }

} // namespace BAM
} // namespace PacBio

#endif // CIGAR_H