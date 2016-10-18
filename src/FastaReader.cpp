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
/// \file FastaReader.cpp
/// \brief Implements the FastaReader class.
//
// Author: Derek Barnett

#include "pbbam/FastaReader.h"
#include <htslib/faidx.h>
#include <fstream>
#include <iostream>
#include <limits>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

struct FastaReaderPrivate
{
    ifstream stream_;
    string name_;
    string bases_;

    FastaReaderPrivate(const std::string& fn)
        : stream_(fn)
    {
        if (!stream_)
            throw std::runtime_error("FastaReader - could not open " + fn + " for reading");
        FetchNext();
    }

    bool GetNext(FastaSequence& record)
    {
        if (name_.empty() && bases_.empty())
            return false;
        record = FastaSequence { name_, bases_ };
        FetchNext();
        return true;
    }

private:
    void FetchNext(void)
    {
        name_.clear();
        bases_.clear();

        SkipNewlines();
        ReadName();
        ReadBases();
    }

    inline void SkipNewlines(void)
    {
        if (!stream_)
            return;
        if (stream_.peek() == '\n')
            stream_.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    void ReadName(void) {
        if (!stream_)
            return;
        if (stream_.get() == '>')
            std::getline(stream_, name_,  '\n');
    }

    void ReadBases(void)
    {
        if (!stream_)
            return;
        char c = static_cast<char>(stream_.peek());
        string line;
        while (c != '>') {
            if (!stream_)
                return;
            std::getline(stream_, line, '\n');
            bases_ += line;
            if (!stream_)
                return;
            c = static_cast<char>(stream_.peek());
        }
    }
};

} // namespace internal
} // namespace BAM
} // namespace PacBio

FastaReader::FastaReader(const std::string& fn)
    : d_{ new internal::FastaReaderPrivate{ fn } }
{ }

FastaReader::FastaReader(FastaReader&& other)
    : d_{ std::move(other.d_) }
{ }

FastaReader& FastaReader::operator=(FastaReader&& other)
{
    d_.swap(other.d_);
    return *this;
}

FastaReader::~FastaReader(void) { }

bool FastaReader::GetNext(FastaSequence& record)
{ return d_->GetNext(record); }

vector<FastaSequence> FastaReader::ReadAll(const string& fn)
{
    vector<FastaSequence> result;
    result.reserve(256);
    FastaReader reader{ fn };
    FastaSequence s;
    while(reader.GetNext(s))
        result.emplace_back(s);
    return result;
}
