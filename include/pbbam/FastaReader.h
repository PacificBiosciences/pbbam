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
/// \file FastaReader.h
/// \brief Defines the FastaReader class.
//
// Author: Derek Barnett

#ifndef FASTAREADER_H
#define FASTAREADER_H

#include "pbbam/FastaSequence.h"
#include <memory>
#include <vector>

namespace PacBio {
namespace BAM {

namespace internal { struct FastaReaderPrivate; }

///
/// \brief The FastaReader provides sequential access to FASTA records.
///
class FastaReader
{
public:
    ///
    /// \brief Reads all FASTA sequences from a file
    ///
    /// \param fn   FASTA filename
    /// \return vector of FastaSequence results
    ///
    static std::vector<FastaSequence> ReadAll(const std::string& fn);

public:
    /// \name Constructors & Related Methods
    /// \{

    explicit FastaReader(const std::string& fn);
    FastaReader(FastaReader&&) = default;
    FastaReader& operator=(FastaReader&&) = default;
    ~FastaReader();

    // copy is disabled
    FastaReader(const FastaReader&) = delete;
    FastaReader& operator=(const FastaReader&) = delete;

    /// \}

public:
    /// \name Sequence Access
    /// \{

    ///
    /// \brief GetNext
    ///
    /// \code{cpp}
    ///
    /// FastaReader reader{ fn };
    /// FastaSequence f;
    /// while (reader.GetNext(f)) {
    ///     // do stuff with f
    /// }
    /// \endcode
    ///
    /// \param[out] record
    /// \return success/failure
    ///
    bool GetNext(FastaSequence& record);

    /// \}

private:
    std::unique_ptr<internal::FastaReaderPrivate> d_;
};

} // namespace BAM
} // namespace PacBio

#endif // FASTAREADER_H
