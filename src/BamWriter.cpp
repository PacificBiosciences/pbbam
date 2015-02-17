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
#include "pbbam/BamRecord.h"
#include "pbbam/BamFile.h"
#include "AssertUtils.h"
#include "MemoryUtils.h"
#include <htslib/bgzf.h>
#include <htslib/hts.h>
#include <iostream>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

BamWriter::BamWriter(const std::string& filename,
                     const SamHeader& header,
                     const BamWriter::CompressionLevel compressionLevel)
    : file_(nullptr)
    , header_(nullptr)
    , error_(BamWriter::NoError)
{
     Open(filename, header.CreateRawData(), compressionLevel);
}

BamWriter::~BamWriter(void)
{
    Close();
}

BamWriter::operator bool(void) const
{
    return (error_ == BamWriter::NoError);
}

void BamWriter::Close(void)
{
    header_.reset();
    file_.reset();
    error_ = BamWriter::NoError;
    filename_.empty();
}

BamWriter::WriteError BamWriter::Error(void) const
{
    return error_;
}

bool BamWriter::Flush(void)
{
    // TODO: sanity checks on file_ & fp
    return (bgzf_flush(file_.get()->fp.bgzf) == 0);
}

bool BamWriter::Open(const string& filename,
                     const shared_ptr<bam_hdr_t> rawHeader,
                     const CompressionLevel compressionLevel)
{
    // ensure clean slate
    Close();

    // store filename
    filename_ = filename;

    // store header
    header_ = rawHeader;
    if (!header_) {
        error_ = BamWriter::NullHeaderError;
        return false;
    }
    PB_ASSERT_OR_RETURN_VALUE((bool)header_, false);

    // open file
    const string& mode = string("wb") + to_string(static_cast<int>(compressionLevel));
    file_.reset(sam_open(filename_.c_str(), mode.c_str()), internal::HtslibFileDeleter());
    if (!file_) {
        error_ = BamWriter::OpenFileError;
        return false;
    }

    // TODO: setup multithreading ??
//    hts_set_threads(file_.get(), 4);

    // write header
    const int ret = sam_hdr_write(file_.get(), header_.get());
    if (ret != 0) {
        error_ = BamWriter::WriteHeaderError;
        return false;
    }

    // if we get here, return success
    return true;
}

bool BamWriter::Write(const BamRecord& record)
{
    return Write(record.RawData());
}

bool BamWriter::Write(const std::shared_ptr<bam1_t>& rawRecord)
{
    const int ret = sam_write1(file_.get(), header_.get(), rawRecord.get());
    if (ret > 0)
        return true;
    else {
        error_ = BamWriter::WriteRecordError;
        return false;
    }
}
