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
/// \file MD5.cpp
/// \brief Implements basic MD5 hash utilities
//
// Author: Brett Bowman

#include "PbbamInternalConfig.h"

#include "pbbam/MD5.h"

#include <stdexcept>

#include <htslib/hts.h>

namespace PacBio {
namespace BAM {

class Md5ContextHelper
{
public:
    Md5ContextHelper() : data_(hts_md5_init())
    {
        if (data_ == nullptr) throw std::runtime_error("could not initialize MD5 context");
    }

    ~Md5ContextHelper() { hts_md5_destroy(data_); }

public:
    std::string Encoded(const std::string& str)
    {
        hts_md5_update(data_, reinterpret_cast<void*>(const_cast<char*>(str.c_str())), str.size());

        unsigned char digest[16];
        hts_md5_final(digest, data_);

        char hexdigest[33];  // leave space for null-term
        hts_md5_hex(hexdigest, digest);

        return std::string{hexdigest, 32};
    }

private:
    hts_md5_context* data_;
};

/// \brief MD5 hash of a string as a 32-digit hexadecimal string
///
std::string MD5Hash(const std::string& str)
{
    Md5ContextHelper md5;
    return md5.Encoded(str);
}

}  // namespace BAM
}  // namespace PacBio
