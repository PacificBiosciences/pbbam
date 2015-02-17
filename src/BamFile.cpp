// Copyright (c) 2014, Pacific Biosciences of California, Inc.
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
#include "AssertUtils.h"
#include "MemoryUtils.h"
#include <htslib/sam.h>
#include <memory>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

BamFile::BamFile(const std::string& filename)
    : filename_(filename)
    , error_(BamFile::NoError)
{
    std::unique_ptr<samFile, internal::HtslibFileDeleter> f(sam_open(filename.c_str(), "rb"));
    if (!f) {
        error_ = BamFile::OpenError;
        return;
    }

    std::shared_ptr<bam_hdr_t> hdr(sam_hdr_read(f.get()), internal::HtslibHeaderDeleter());
    if (!hdr) {
        error_ = BamFile::ReadHeaderError;
        return;
    }
    header_ = std::move(SamHeader::FromRawData(hdr));
}

BamFile::BamFile(const BamFile& other)
    : filename_(other.filename_)
    , header_(other.header_)
    , error_(other.error_)
{ }

BamFile::~BamFile(void) { }

BamFile::operator bool(void) const
{
    return error_ == BamFile::NoError;
}

BamFile::FileError BamFile::Error(void) const
{
    return error_;
}

string BamFile::Filename(void) const
{
    return filename_;
}

SamHeader BamFile::Header(void) const
{
    return header_;
}

bool BamFile::IsPacBioBAM(void) const
{
    return !header_.pacbioBamVersion.empty();
}

string BamFile::StandardIndexFilename(void) const
{
    return filename_ + ".bai";
}

string BamFile::PacBioIndexFilename(void) const
{
    return filename_ + ".pbi";
}

int BamFile::ReferenceId(const string& name) const
{
    return header_.sequences.IndexOf(name);
}

uint32_t BamFile::ReferenceLength(const std::string& name) const
{
    return ReferenceLength(ReferenceId(name));
}

uint32_t BamFile::ReferenceLength(const int id) const
{
    if (id < 0)
        return 0;
    try {
        return stoul(header_.sequences.At(id).length);
    } catch (std::exception&) {
        return 0;
    }
}

string BamFile::ReferenceName(const int id) const
{
    if (id < 0)
        return string();
    try {
        return header_.sequences.At(id).name;
    } catch (std::exception&) {
        return string();
    }
}
