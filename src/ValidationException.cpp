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
/// \file ValidationException.cpp
/// \brief Implements the ValidationException class.
//
// Author: Derek Barnett

#include "pbbam/exception/ValidationException.h"
using namespace PacBio;
using namespace PacBio::BAM;

ValidationException::ValidationException(const ErrorMap& fileErrors,
                                         const ErrorMap& readGroupErrors,
                                         const ErrorMap& recordErrors)
    : std::runtime_error("")
    , fileErrors_(fileErrors)
    , readGroupErrors_(readGroupErrors)
    , recordErrors_(recordErrors)
{
    FormatMessage();
}

ValidationException::ValidationException(ErrorMap&& fileErrors,
                                         ErrorMap&& readGroupErrors,
                                         ErrorMap&& recordErrors)
    : std::runtime_error("")
    , fileErrors_(std::move(fileErrors))
    , readGroupErrors_(std::move(readGroupErrors))
    , recordErrors_(std::move(recordErrors))
{
    FormatMessage();
}

const ValidationException::ErrorMap& ValidationException::FileErrors(void) const
{ return fileErrors_; }

const ValidationException::ErrorMap& ValidationException::ReadGroupErrors(void) const
{ return readGroupErrors_; }

const ValidationException::ErrorMap& ValidationException::RecordErrors(void) const
{ return recordErrors_; }

const char* ValidationException::what(void) const noexcept
{ return msg_.c_str(); }

void ValidationException::FormatMessage(void)
{
    std::stringstream s;
    s << "Validation failed: " << std::endl;

    // file errors
    if (!fileErrors_.empty()) {
        auto fileIter = fileErrors_.cbegin();
        auto fileEnd  = fileErrors_.cend();
        for ( ; fileIter != fileEnd; ++fileIter) {
            s << "  In file (" << fileIter->first << ") : " << std::endl;
            const auto& errors = fileIter->second;
            for (const auto& e : errors)
                s << "    " << e << std::endl;
        }
    }

    // read group errors
    if (!readGroupErrors_.empty()) {
        auto rgIter = readGroupErrors_.cbegin();
        auto rgEnd  = readGroupErrors_.cend();
        for ( ; rgIter != rgEnd; ++rgIter) {
            s << "  In read group (" << rgIter->first << ") : " << std::endl;
            const auto& errors = rgIter->second;
            for (const auto& e : errors)
                s << "    " << e << std::endl;
        }
    }

    // record errors
    if (!recordErrors_.empty()) {
        auto recIter = recordErrors_.cbegin();
        auto recEnd  = recordErrors_.cend();
        for ( ; recIter != recEnd; ++recIter) {
            s << "  In record (" << recIter->first << ") : " << std::endl;
            const auto& errors = recIter->second;
            for (const auto& e : errors)
                s << "    " << e << std::endl;
        }
    }

    msg_ = s.str();
}
