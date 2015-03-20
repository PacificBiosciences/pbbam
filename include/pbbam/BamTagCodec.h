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

#ifndef BAMTAGCODEC_H
#define BAMTAGCODEC_H

#include "pbbam/Config.h"
#include "pbbam/TagCollection.h"
#include <vector>

namespace PacBio {
namespace BAM {

class PBBAM_EXPORT BamTagCodec
{

// high-level, operate on a full collection
public:
    static TagCollection Decode(const std::vector<uint8_t>& data);
    static std::vector<uint8_t> Encode(const PacBio::BAM::TagCollection& tags);

// per-tag methods
public:

    // returns the SAM/BAM single char code for tag type
    static uint8_t TagTypeCode(const PacBio::BAM::Tag& tag);

    // returns the tag value's raw data in bytes
    // NOTE: does *NOT* encode name & tag type. It does however,
    // include the element type of an array tag
    static std::vector<uint8_t> ToRawData(const PacBio::BAM::Tag& tag);

    // TODO: make this hidden a bit more, maybe this whole class in fact
    // rawData should be the result of sam.h:bam_aux_get(...)
    static PacBio::BAM::Tag FromRawData(uint8_t* rawData);
};

} // namespace BAM
} // namespace PacBio

#endif // BAMTAGCODEC_H
