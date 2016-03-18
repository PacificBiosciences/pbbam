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
/// \file ValidationErrors.cpp
/// \brief Implements the ValidationErrors class.
//
// Author: Derek Barnett

#include "ValidationErrors.h"
#include "pbbam/exception/ValidationException.h"
#include "StringUtils.h"
#include <sstream>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace std;

const size_t ValidationErrors::MAX;

ValidationErrors::ValidationErrors(const size_t maxNumErrors)
    : maxNumErrors_(maxNumErrors)
    , currentNumErrors_(0)
{
    if (maxNumErrors_ == 0)
        maxNumErrors_ = ValidationErrors::MAX;
}

void ValidationErrors::AddFileError(const std::string& fn,
                                    const std::string& details)
{
    string copy = details;
    AddFileError(fn, std::move(copy));
}

void ValidationErrors::AddFileError(const std::string& fn,
                                    std::string&& details)
{
    fileErrors_[fn].push_back(std::move(details));
    OnErrorAdded();
}

void ValidationErrors::AddReadGroupError(const std::string& rg,
                                         const std::string& details)
{
    string copy = details;
    AddReadGroupError(rg, std::move(copy));
}

void ValidationErrors::AddReadGroupError(const std::string& rg,
                                         std::string&& details)
{
    readGroupErrors_[rg].push_back(std::move(details));
    OnErrorAdded();
}

void ValidationErrors::AddRecordError(const std::string& name,
                                      const std::string& details)
{
    string copy = details;
    AddRecordError(name, std::move(copy));
}

void ValidationErrors::AddRecordError(const std::string& name,
                                      std::string&& details)
{
    recordErrors_[name].push_back(std::move(details));
    OnErrorAdded();
}

void ValidationErrors::AddTagLengthError(const string& name,
                                         const string& tagLabel,
                                         const string& tagName,
                                         const size_t observed,
                                         const size_t expected)
{
    string copy  = tagLabel;
    string copy2 = tagName;
    AddTagLengthError(name, std::move(copy), std::move(copy2), observed, expected);
}

void ValidationErrors::AddTagLengthError(const string& name,
                                         string&& tagLabel,
                                         string&& tagName,
                                         const size_t observed,
                                         const size_t expected)
{
    // format
    stringstream s;
    s << tagLabel << " tag (" << tagName << ") length: " << observed
      << ", does not match expected length: " << expected;
    AddRecordError(name, s.str());
}

bool ValidationErrors::IsEmpty(void) const
{
    return currentNumErrors_ == 0;
}

void ValidationErrors::OnErrorAdded(void)
{
    ++currentNumErrors_;
    if (currentNumErrors_ == maxNumErrors_)
        ThrowErrors();
}

void ValidationErrors::ThrowErrors(void)
{
    throw ValidationException(std::move(fileErrors_),
                              std::move(readGroupErrors_),
                              std::move(recordErrors_));
}
