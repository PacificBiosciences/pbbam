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

#ifndef SAMHEADER_H
#define SAMHEADER_H

#include "htslib/sam.h"
#include "pbbam/Config.h"
#include "pbbam/DictionaryBase.h"
#include "pbbam/SamProgram.h"
#include "pbbam/SamReadGroup.h"
#include "pbbam/SamSequence.h"
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

class BamReader;

class PBBAM_EXPORT SamHeader
{
public:
    typedef DictionaryBase<SamReadGroup> ReadGroupDictionary;
    typedef DictionaryBase<SamSequence>  SequenceDictionary;
    typedef DictionaryBase<SamProgram>   ProgramDictionary;

public:
    SamHeader(void);
    SamHeader(const SamHeader& other) = default;
    SamHeader(SamHeader&& other) = default;
    ~SamHeader(void) { }

public:
    std::string version;
    std::string sortOrder;
    std::string pacbioBamVersion;
    ReadGroupDictionary readGroups;
    SequenceDictionary  sequences;
    ProgramDictionary   programs;
    std::vector<std::string> comments;

#ifdef PBBAM_TESTING
public:
#else
private:
#endif // PBBAM_TESTING

    // converts SamHeader contents to htslib raw data
    // NOTE: caller takes ownership of object (so it must be cleaned later up via sam.h:bam_hdr_destroy())
    bam_hdr_t* CreateRawData(void) const;

    // returns a SamHeader object, with a deep copy of @rawData contents
    static SamHeader FromRawData(bam_hdr_t* rawData);

private:
    friend class BamReader;
    friend class BamWriter;
};

} // namespace BAM
} // namespace PacBio

#endif // SAMHEADER_H
