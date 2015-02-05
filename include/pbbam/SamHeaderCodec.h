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

#ifndef SAMHEADERCODEC_H
#define SAMHEADERCODEC_H

#include "pbbam/Config.h"
#include "pbbam/SamHeader.h"
#include <sstream>
#include <string>

namespace PacBio {
namespace BAM {

class PBBAM_EXPORT SamHeaderCodec
{
public:
    static SamHeader Decode(const std::string& text);
    static std::string Encode(const SamHeader& header);

private:

    // decode helper functions
    static void DecodeHeaderLine(const std::string& line, SamHeader* header);
    static void DecodeSequenceLine(const std::string& line, SamHeader* header);
    static void DecodeReadGroupLine(const std::string& line, SamHeader* header);
    static void DecodeProgramLine(const std::string& line, SamHeader* header);
    static void DecodeCommentLine(const std::string& line, SamHeader* header);

    // encode helper functions
    static void EncodeHeaderLine(const std::string& version,
                                 const std::string& sortOrder,
                                 const std::string& pbVersion,
                                 std::stringstream* out);
    static void EncodeSequence(const SamSequence& sequence, std::stringstream* out);
    static void EncodeReadGroup(const SamReadGroup& readGroup, std::stringstream* out);
    static void EncodeProgram(const SamProgram& program, std::stringstream* out);
    static void EncodeComment(const std::string& comment, std::stringstream* out);
};

} // namespace BAM
} // namespace PacBio

#endif // SAMHEADERCODEC_H
