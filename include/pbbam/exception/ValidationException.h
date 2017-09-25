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
/// \file ValidationException.h
/// \brief Defines the ValidationException class.
//
// Author: Derek Barnett

#ifndef VALIDATIONEXCEPTION_H
#define VALIDATIONEXCEPTION_H

#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

/// \brief The ValidationExecption represents an exception that will be thrown
///        when any error is encountered using the Validator API. In addition to
///        a default display message, it provides programmatic access to all
///        reported error messages.
///
/// \sa Validator::Validate(const BamRecord& record)
///
class ValidationException : public std::runtime_error
{
public:
    using ErrorList = std::vector<std::string>;
    using ErrorMap  = std::map<std::string, ErrorList>;

public:
    ValidationException(ErrorMap fileErrors,
                        ErrorMap readGroupErrors,
                        ErrorMap recordErrors);

    // This is a work around for the Intel PHI compiler (icpc)
    ~ValidationException() throw()
    {
    }

public:
    const ErrorMap& FileErrors() const;
    const ErrorMap& ReadGroupErrors() const;
    const ErrorMap& RecordErrors() const;

    const char* what() const noexcept override;

private:
    ErrorMap fileErrors_;
    ErrorMap readGroupErrors_;
    ErrorMap recordErrors_;
    std::string msg_;

private:
    void FormatMessage();
};

} // namespace BAM
} // namespace PacBio

#endif // VALIDATIONEXCEPTION_H
