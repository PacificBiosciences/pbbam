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
/// \file InvalidSequencingChemistryException.h
/// \brief Defines the InvalidSequencingChemistryException class.
//
// Author: Derek Barnett

#ifndef INVALIDSEQUENCINGCHEMISTRYEXCEPTION_H
#define INVALIDSEQUENCINGCHEMISTRYEXCEPTION_H

#include <exception>
#include <sstream>
#include <string>

namespace PacBio {
namespace BAM {

/// \brief The InvalidSequencingChemistryException class represents an exception
///        that will be thrown when an invalid sequencing chemistry combination
///        is encountered.
///
class InvalidSequencingChemistryException : public std::exception
{
public:
    InvalidSequencingChemistryException(std::string bindingKit,
                                        std::string sequencingKit,
                                        std::string basecallerVersion)
        : bindingKit_(std::move(bindingKit))
        , sequencingKit_(std::move(sequencingKit))
        , basecallerVersion_(std::move(basecallerVersion))
    {
        std::stringstream s;
        s << "unsupported sequencing chemistry combination: " << std::endl
          << "    binding kit:        " << bindingKit_ << std::endl
          << "    sequencing kit:     " << sequencingKit_ << std::endl
          << "    basecaller version: " << basecallerVersion_ << std::endl;
        what_ = s.str();
    }

    // This is a work around for the Intel PHI compiler (icpc)
    ~InvalidSequencingChemistryException() throw()
    {
    }

public:
    const std::string& BindingKit() const
    { return bindingKit_; }

    const std::string& SequencingKit() const
    { return sequencingKit_; }

    const std::string& BasecallerVersion() const
    { return basecallerVersion_; }

public:
    const char* what() const noexcept override
    { return what_.c_str(); }

protected:
    std::string bindingKit_;
    std::string sequencingKit_;
    std::string basecallerVersion_;
    std::string what_;
};

} // namespace BAM
} // namespace PacBio

#endif // INVALIDSEQUENCINGCHEMISTRYEXCEPTION_H
