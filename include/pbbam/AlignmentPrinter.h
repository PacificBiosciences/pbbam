// Copyright (c) 2015, Pacific Biosciences of California, Inc.
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
/// \file AlignmentPrinter.h
/// \brief Defines the AlignmentPrinter class.
//
// Author: Armin Töpfer

#ifndef ALIGNMENTPRINTER_H
#define ALIGNMENTPRINTER_H

#include <memory>
#include <string>
#include "pbbam/BamRecord.h"
#include "pbbam/IndexedFastaReader.h"
#include "pbbam/Orientation.h"

namespace PacBio {
namespace BAM {

class BamRecord;

/// \brief The AlignmentPrinter class "pretty-prints" an alignment with respect
///        to its associated reference sequence.
///
/// Example output:
/// \verbinclude plaintext/AlignmentPrinterOutput.txt
///
class AlignmentPrinter
{
public:
    /// \name Constructors & Related Methods
    /// \{

    /// Constructs the alignment printer with an associated FASTA file reader.
    ///
    /// \param[in] ifr FASTA reader
    ///
    /// \throws std::runtime_error if FASTA file cannot be opened for reading.
    ///
    AlignmentPrinter(const IndexedFastaReader& ifr);

    AlignmentPrinter() = delete;
    AlignmentPrinter(const AlignmentPrinter&) = delete;
    AlignmentPrinter(AlignmentPrinter&&) = default;
    AlignmentPrinter& operator=(const AlignmentPrinter&) = delete;
    AlignmentPrinter& operator=(AlignmentPrinter&&) = default;
    ~AlignmentPrinter() = default;

    /// \}

public:
    /// \name Printing
    /// \{

    /// Pretty-prints an aligned BamRecord to std::string.
    ///
    /// \note The current implementation includes ANSI escape sequences for
    ///       coloring terminal output. Future versions of this method will
    ///       likely make this optional.
    ///
    /// \returns formatted string containing the alignment and summary
    ///          information
    ///
    std::string Print(const BamRecord& record,
                      const Orientation orientation = Orientation::GENOMIC);

    /// \}

private:
	const std::unique_ptr<IndexedFastaReader> ifr_;
};

} // namespace BAM
} // namespace PacBio

#endif // ALIGNMENTPRINTER_H
