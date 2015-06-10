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

// Author: David Alexander

#ifndef INDEXEDFASTAREADER_H
#define INDEXEDFASTAREADER_H

#include "pbbam/Orientation.h"
#include "pbbam/Position.h"
#include "htslib/faidx.h"

#include <string>
#include <iostream>
#include <stdexcept>

namespace PacBio {
namespace BAM {

class GenomicInterval;
class BamRecord;

class IndexedFastaReader {

public:
    IndexedFastaReader() = delete;
    IndexedFastaReader(const std::string& filename);
    ~IndexedFastaReader();

public:
    // Copy constructor
    IndexedFastaReader(const IndexedFastaReader& src)
    {
        if (!Open(src.filename_)) 
            throw std::runtime_error("Cannot open file " + src.filename_);
    }

    // Copy assignment operator
    IndexedFastaReader& operator=(const IndexedFastaReader& rhs)
    {
        if(&rhs == this) return *this;

        Open(rhs.filename_);
        return *this;
    }

public:
    std::string Subsequence(const std::string& id, Position begin, Position end) const;
    std::string Subsequence(const GenomicInterval& interval) const;
    std::string Subsequence(const char* htslibRegion) const;

public:
    // \returns subsequence of the reference corresponding to the bamRecord,
    // oriented and gapped as requested.  For example, "native" orientation
    // and "gapped" will return the reference sequence with gaps inserted, as
    // would align against the read in "native" orientation
    std::string ReferenceSubsequence(const BamRecord& bamRecord,
                                     const Orientation orientation=Orientation::GENOMIC,
                                     const bool gapped=false,
                                     const bool exciseSoftClips=false) const;

public:
    int NumSequences() const;
    bool HasSequence(const std::string& name) const;
    int SequenceLength(const std::string& name) const;

private:
    std::string filename_;
    faidx_t* handle_;

private:
    void Close(void);
    bool Open(const std::string& filename);
};



}  // PacBio
}  // BAM
#endif  // INDEXEDFASTAREADER_H
