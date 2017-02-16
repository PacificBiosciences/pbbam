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
/// \file CompositeFastaReader.h
/// \brief Defines the composite FASTA reader, for working with multiple input
///       files.
//
// Author: Derek Barnett

#ifndef COMPOSITEFASTAREADER_H
#define COMPOSITEFASTAREADER_H

#include "pbbam/Config.h"
#include "pbbam/DataSet.h"
#include "pbbam/FastaReader.h"
#include "pbbam/FastaSequence.h"

#include <deque>
#include <memory>
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

/// \brief The CompositeFastaReader class provides read access to
///        multiple FASTA files, reading through the entire contents of each
///        file.
///
/// Input files will be accessed in the order provided to the constructor. Each
/// file's contents will be exhausted before moving on to the next one (as
/// opposed to a "round-robin" scheme).
///
class PBBAM_EXPORT CompositeFastaReader
{
public:
    /// \name Contstructors & Related Methods
    /// \{

    CompositeFastaReader(const std::vector<std::string>& fastaFiles);
    CompositeFastaReader(std::vector<std::string>&& fastaFiles);
    CompositeFastaReader(const DataSet& dataset);

    /// \}

public:
    /// \name Data Access
    /// \{

    /// Fetches next FASTA sequence.
    ///
    /// \returns true on success, false if no more data available.
    ///
    bool GetNext(FastaSequence& seq);

    /// \}

private:
    std::deque<std::unique_ptr<FastaReader> > readers_;
};

} // namespace BAM
} // namespace PacBio

#include "internal/CompositeFastaReader.inl"

#endif // COMPOSITEFASTAREADER_H
