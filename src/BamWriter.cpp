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

#include "pbbam/BamWriter.h"
#include "AssertUtils.h"
#include "MemoryUtils.h"
#include <htslib/hts.h>
#include <iostream>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

BamWriter::BamWriter()
    : file_(nullptr)
    , header_(nullptr)
{ }

BamWriter::~BamWriter(void)
{
    Close();
}

void BamWriter::Close(void)
{
    filename_.clear();
    errorString_.clear();

    file_.reset();
    header_.reset();
}

string BamWriter::ErrorString(void) const
{
    return errorString_;
}

bool BamWriter::HasError(void) const
{
    return !errorString_.empty();
}

bool BamWriter::Open(const string& filename,
                     bam_hdr_t* rawHeader,
                     const CompressionLevel compressionLevel)
{
    // ensure clean slate
    Close();

    // store filename
    filename_ = filename;

    // store header
    header_.reset(rawHeader, internal::RawHeaderDeleter());
    if (!header_) {
        errorString_ = "invalid header data";
        return false;
    }

    PB_ASSERT_OR_RETURN_VALUE(header_, false);

    // open file
    const string& mode = string("wb") + to_string(static_cast<int>(compressionLevel));
    file_.reset(sam_open(filename_.c_str(), mode.c_str()), internal::RawFileDeleter());
    if (!file_) {
        errorString_ = "could not open " + filename_;
        return false;
    }

//    hts_set_threads(file_.get(), 4);

    // write header
    const int ret = sam_hdr_write(file_.get(), header_.get());
    if (ret != 0) {
        errorString_ = "could not write header";
        return false;
    }

    return true;
}

bool BamWriter::Open(const string& filename,
                     const SamHeader& header,
                     const CompressionLevel compressionLevel)
{
    return Open(filename, header.CreateRawData(), compressionLevel);
}

bool BamWriter::Write(const BamRecord& record)
{
    return Write(record.RawData().get());
}

bool BamWriter::Write(const bam1_t* rawRecord)
{
    const int ret = sam_write1(file_.get(), header_.get(), rawRecord);
    if (ret > 0)
        return true;
    else {
        errorString_ = "could not write BAM record";
        return false;
    }
}
