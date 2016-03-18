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
/// \file ValidationErrors.h
/// \brief Defines the ValidationErrors class.
//
// Author: Derek Barnett

#ifndef VALIDATIONERRORS_H
#define VALIDATIONERRORS_H

#include <limits>
#include <map>
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {
namespace internal {

/// The ValidationErrors class catches error messages accumulated during
/// validation (see Validator).
///
/// Convenience methods are provided for different BAM components, to help
/// format the displayed output.
///
/// A maximum number of errors can be provided at construction, and this class
/// will automatially throw a ValidationException whenever that count is reached.
/// Otherwise, the Validator will check IsEmpty() and call ThrowErrors() if true.
///
class ValidationErrors
{
public:
    typedef std::vector<std::string>         ErrorList;
    typedef std::map<std::string, ErrorList> ErrorMap;
public:
    static const size_t MAX = std::numeric_limits<size_t>::max();

public:
    ValidationErrors(const size_t maxNumErrors = ValidationErrors::MAX);

public:
    void AddFileError(const std::string& fn, const std::string& details);
    void AddFileError(const std::string& fn, std::string&& details);

    void AddReadGroupError(const std::string& rg, const std::string& details);
    void AddReadGroupError(const std::string& rg, std::string&& details);

    void AddRecordError(const std::string& name, const std::string& details);
    void AddRecordError(const std::string& name, std::string&& details);

    void AddTagLengthError(const std::string& name,
                           const std::string& tagLabel,
                           const std::string& tagName,
                           const size_t observed,
                           const size_t expected);
    void AddTagLengthError(const std::string& name,
                           std::string&& tagLabel,
                           std::string&& tagName,
                           const size_t observed,
                           const size_t expected);

public:
    bool IsEmpty(void) const;
    void ThrowErrors(void);

private:
    size_t maxNumErrors_;
    size_t currentNumErrors_;
    ErrorMap fileErrors_;
    ErrorMap readGroupErrors_;
    ErrorMap recordErrors_;

private:
    void OnErrorAdded(void);
};

} // namespace internal
} // namespace BAM
} // namespace PacBio

#endif // VALIDATIONERRORS_H
