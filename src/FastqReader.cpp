// Copyright (c) 2016, Pacific Biosciences of California, Inc.
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
/// \file FastqReader.cpp
/// \brief Implements the FastqReader class.
//
// Author: Derek Barnett

#include "pbbam/FastqReader.h"
#include <htslib/faidx.h>
#include <stdexcept>
#include <fstream>
#include <iostream>
#include <limits>

namespace PacBio {
namespace BAM {
namespace internal {

struct FastqReaderPrivate {
public:
    explicit FastqReaderPrivate(const std::string& fn)
        : stream_(fn)
    {
        if (!stream_)
            throw std::runtime_error("FastqReader - could not open " + fn + " for reading");
        FetchNext();
    }

    bool GetNext(FastqSequence& record)
    {
        if (name_.empty() && bases_.empty() && quals_.empty())
            return false;
        record = FastqSequence { name_, bases_, quals_ };
        FetchNext();
        return true;
    }

private:
    void FetchNext(void)
    {
        name_.clear();
        bases_.clear();
        quals_.clear();

        if (!stream_ || stream_.eof())
            return;

        SkipNewlines();

        ReadName();
        ReadBases();
        stream_.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // ignore "comment line"
        ReadQuals();
    }

    inline void SkipNewlines(void)
     {
         if (stream_.peek() == '\n')
             stream_.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
     }

     void ReadName(void) {
         if (stream_.get() == '@')
             std::getline(stream_, name_,  '\n');
     }

     void ReadBases(void)
     {
         std::getline(stream_, bases_,  '\n');
     }

     void ReadQuals(void)
     {
         std::getline(stream_, quals_,  '\n');
     }

 private:
     std::ifstream stream_;
     std::string name_;
     std::string bases_;
     std::string quals_;
};

} // namespace internal

FastqReader::FastqReader(const std::string& fn)
    : d_{ new internal::FastqReaderPrivate{ fn } }
{ }

FastqReader::FastqReader(FastqReader&& other)
    : d_{ std::move(other.d_) }
{ }

FastqReader& FastqReader::operator=(FastqReader&& other)
{
    if (this != &other) {
        d_.swap(other.d_);
    }
    return *this;
}

FastqReader::~FastqReader(void) { }

bool FastqReader::GetNext(FastqSequence& record)
{ return d_->GetNext(record); }

std::vector<FastqSequence> FastqReader::ReadAll(const std::string& fn)
{
    std::vector<FastqSequence> result;
    result.reserve(256);
    FastqReader reader{ fn };
    FastqSequence s;
    while(reader.GetNext(s))
        result.emplace_back(s);
    return result;
}

} // namespace BAM
} // namespace PacBio
