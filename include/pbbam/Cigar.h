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
/// \file Cigar.h
/// \brief Defines the Cigar class.
//
// Author: Derek Barnett

#ifndef CIGAR_H
#define CIGAR_H

#include "pbbam/CigarOperation.h"
#include "pbbam/Config.h"
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

/// \brief The Cigar class represents the CIGAR string used to report alignment
///        charateristics in SAM/BAM.
///
/// \note Use of the 'M' operator is forbidden in PacBio BAMs. See
///       CigarOperationType description for more information.
///
/// \sa https://samtools.github.io/hts-specs/SAMv1.pdf for more information on CIGAR in general.
///
class PBBAM_EXPORT Cigar : public std::vector<CigarOperation>
{
public:
    /// \name Constructors & Related Methods
    /// \{

    /// \brief Creates a Cigar object from SAM/BAM string input
    ///
    /// \param [in] stdString   SAM/BAM formatted CIGAR data
    /// \returns a Cigar object representing the input data
    ///
    /// \note This class may be removed from the public API in the future,
    ///       as the constructor taking a std::string accomplishes the same end.
    ///
    static Cigar FromStdString(const std::string& stdString);

    /// \brief Creates an empty Cigar.
    Cigar() = default;

    /// \brief Creates a Cigar object from SAM/BAM string input
    ///
    /// \param [in] cigarString   SAM/BAM formatted CIGAR data
    ///
    Cigar(const std::string& cigarString);

    Cigar(const Cigar& other) = default;
    Cigar(Cigar&& other) = default;
    Cigar& operator=(const Cigar& other) = default;
    Cigar& operator=(Cigar&& other) = default;
    ~Cigar() = default;

    /// \}

public:
    /// \name Conversion Methods
    /// \{

    /// Converts Cigar object data to SAM/BAM formatted string
    ///
    /// \returns SAM/BAM formatted std::string
    ///
    std::string ToStdString() const;

    /// \}
};

} // namespace BAM
} // namespace PacBio

#include "pbbam/internal/Cigar.inl"

#endif // CIGAR_H
