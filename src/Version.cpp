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
/// \file Version.cpp
/// \brief Implements the Version class.
//
// Author: Derek Barnett

#include "Version.h"
#include "SequenceUtils.h"
#include <sstream>

using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace std;

const Version Version::Current = Version(3,0,3);
const Version Version::Minimum = Version(3,0,1);

// string must be "<major>.<minor>.<version>"
Version::Version(const std::string& v)
    : major_(0)
    , minor_(0)
    , revision_(0)
{
    // parse string
    try {
        const auto fields = internal::Split(v, '.');
        const auto numFields = fields.size();
        if (numFields == 0)
            throw std::runtime_error("invalid version number - empty string");
        major_ = std::stoi(fields.at(0));
        if (numFields > 1) {
            minor_ = std::stoi(fields.at(1));
            if (numFields > 2 )
                revision_ = std::stoi(fields.at(2));
        }
    } catch (std::exception&) {
        auto msg = string{ "invalid version number (" + v + "): failed to parse" };
        throw std::runtime_error(msg);
    }

    // ensure valid numbers
    Check();
}

std::string Version::ToString(void) const
{
    std::stringstream s;
    s << major_ << '.' << minor_ << '.' << revision_;
    return s.str();
}

