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

// Author: Derek Barnett

#include "pbbam/BamFile.h"
#include "MemoryUtils.h"
#include <htslib/sam.h>
#include <memory>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

BamFile::BamFile(void)
    : error_(BamFile::NoError)
    , header_(nullptr)
{ }

BamFile::BamFile(const std::string& filename)
    : error_(BamFile::NoError)
    , header_(nullptr)
{
    Open(filename);
}

BamFile::BamFile(const BamFile& other)
    : filename_(other.filename_)
    , error_(other.error_)
    , header_(other.header_)
{ }

BamFile::~BamFile(void) { }

void BamFile::Close(void)
{
    error_ = BamFile::NoError;
    filename_.clear();
    header_.reset();
}

bool BamFile::IsOpen(void) const
{ return !filename_.empty() && header_; }

void BamFile::Open(const string& filename)
{
    filename_ = filename;

    // attempt open
    std::unique_ptr<samFile, internal::HtslibFileDeleter> f(sam_open(filename.c_str(), "rb"));
    if (!f) {
        error_ = BamFile::OpenError;
        return;
    }

    // attempt fetch header
    std::unique_ptr<bam_hdr_t, internal::HtslibHeaderDeleter> hdr(sam_hdr_read(f.get()));
    if (!hdr) {
        error_ = BamFile::ReadHeaderError;
        return;
    }
    header_ = internal::BamHeaderMemory::FromRawData(hdr.get());
}
