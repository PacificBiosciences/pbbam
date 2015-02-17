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

#include "pbbam/BamReader.h"
#include "MemoryUtils.h"
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

BamReader::BamReader(void)
    : file_(nullptr)
    , header_(nullptr)
    , error_(BamReader::NoError)
{ }

BamReader::~BamReader(void)
{
    Close();
}

void BamReader::Close(void)
{
    filename_.clear();

    file_.reset();
    header_.reset();

    error_ = BamReader::NoError;
}

BamReader::ReadError BamReader::Error(void) const
{
    return error_;
}

bool BamReader::HasError(void) const
{
    return error_ != BamReader::NoError;
}

bool BamReader::GetNext(std::shared_ptr<BamRecord> record)
{
    return GetNext(record->d_);
}

bool BamReader::GetNext(std::shared_ptr<bam1_t> rawRecord)
{
    const int ret = sam_read1(file_.get(), header_.get(), rawRecord.get());
    if (ret >= 0)
        return true;

    // else determine error type
    // -1 : EOF - normal (so should we set error then? or just return false so we're done)
    // -2 : EOF - truncated
    // -3 : malformatted fixed-length data
    // -4 : malformatted variable-length data
    error_ = BamReader::ReadRecordError;
    return false;
}

SamHeader BamReader::Header(void) const {
    return SamHeader::FromRawData(header_);
}

bool BamReader::Open(const string& filename)
{
    // ensure clean slate
    Close();

    // open file
    filename_ = filename;
    file_.reset(sam_open(filename_.c_str(), "rb"), internal::HtslibFileDeleter());
    if (!file_) {
        error_ = BamReader::OpenFileError;
        return false;
    }

    // read header
    header_.reset(sam_hdr_read(file_.get()), internal::HtslibHeaderDeleter());
    if (!header_) {
        error_ = BamReader::ReadHeaderError;
        return false;
    }

    // if we get here, return success
    return true;
}

string BamReader::PacBioBamVersion(void) const {
    const SamHeader& header = SamHeader::FromRawData(header_);
    return header.pacbioBamVersion;
}

std::shared_ptr<bam_hdr_t> BamReader::RawHeader(void) const {
    return header_;
}
