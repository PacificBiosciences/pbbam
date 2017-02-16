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
/// \file CompositeFastaReader.inl
/// \brief Inline implementation for the composite FASTA reader, for
///        working with multiple input files.
//
// Author: Derek Barnett

#include "pbbam/CompositeFastaReader.h"

namespace PacBio {
namespace BAM {


inline CompositeFastaReader::CompositeFastaReader(const std::vector<std::string>& fastaFiles)
{
    for (auto&& fn : fastaFiles)
        readers_.emplace_back(new FastaReader{ fn });
}

inline CompositeFastaReader::CompositeFastaReader(std::vector<std::string>&& fastaFiles)
{
    for (auto&& fn : fastaFiles)
        readers_.emplace_back(new FastaReader{ std::move(fn) });
}

inline CompositeFastaReader::CompositeFastaReader(const DataSet& dataset)
    : CompositeFastaReader(dataset.FastaFiles())
{ }

inline bool CompositeFastaReader::GetNext(FastaSequence& seq)
{
    // try first reader, if successful return true
    // else pop reader and try next, until all readers exhausted
    while (!readers_.empty()) {
        auto& reader = readers_.front();
        if (reader->GetNext(seq))
            return true;
        else
            readers_.pop_front();
    }

    // no readers available
    return false;
}

} // namespace BAM
} // namespace PacBio
