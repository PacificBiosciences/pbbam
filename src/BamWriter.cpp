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

#include "pbbam/BamWriter.h"
#include "pbbam/BamFile.h"
#include "AssertUtils.h"
#include "MemoryUtils.h"
#include <htslib/bgzf.h>
#include <htslib/hts.h>
#include <thread>
#include <iostream>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

class BamWriterPrivate
{
public:
    BamWriterPrivate(void)
        : file_(nullptr)
        , header_(nullptr)
    { }

public:
    void Open(const std::string& filename,
              const PBBAM_SHARED_PTR<bam_hdr_t> rawHeader,
              const BamWriter::CompressionLevel compressionLevel = BamWriter::DefaultCompression,
              size_t numThreads = 4);
    void Write(const PBBAM_SHARED_PTR<bam1_t>& rawRecord);

public:
    std::unique_ptr<samFile, internal::HtslibFileDeleter> file_;
    PBBAM_SHARED_PTR<bam_hdr_t> header_;
    std::string filename_;
};

void BamWriterPrivate::Open(const string& filename,
                            const PBBAM_SHARED_PTR<bam_hdr_t> rawHeader,
                            const BamWriter::CompressionLevel compressionLevel,
                            size_t numThreads)
{
    // store filename
    filename_ = filename;

    // store header
    header_ = rawHeader;
    if (!header_)
        throw std::exception();

    // open file
    const string& mode = string("wb") + to_string(static_cast<int>(compressionLevel));
    file_.reset(sam_open(filename_.c_str(), mode.c_str()));
    if (!file_)
        throw std::exception();

    // if no explicit thread count given, attempt built-in check
    if (numThreads == 0) {
        numThreads = thread::hardware_concurrency();

        // if still unknown, default to single-threaded
        if (numThreads == 0)
            numThreads = 1;
    }

    // if multithreading requested, enable it
    if (numThreads > 1)
        hts_set_threads(file_.get(), numThreads);

    // write header
    const int ret = sam_hdr_write(file_.get(), header_.get());
    if (ret != 0)
        throw std::exception();
}

void BamWriterPrivate::Write(const PBBAM_SHARED_PTR<bam1_t>& rawRecord)
{
    const int ret = sam_write1(file_.get(), header_.get(), rawRecord.get());
    if (ret <= 0)
        throw std::exception();
}

} // namespace internal
} // namespace BAM
} // namespace PacBio

BamWriter::BamWriter(const std::string& filename,
                     const BamHeader& header,
                     const BamWriter::CompressionLevel compressionLevel,
                     const size_t numThreads)
    : d_(new internal::BamWriterPrivate)
{
     d_->Open(filename,
              internal::BamHeaderMemory::MakeRawHeader(header),
              compressionLevel,
              numThreads);
}

BamWriter::~BamWriter(void)
{
    bgzf_flush(d_->file_.get()->fp.bgzf);
}

void BamWriter::TryFlush(void)
{
    // TODO: sanity checks on file_ & fp
    const int ret = bgzf_flush(d_->file_.get()->fp.bgzf);
    if (ret != 0)
        throw std::exception();
}

void BamWriter::Write(const BamRecord& record)
{ d_->Write(internal::BamRecordMemory::GetRawData(record)); }

void BamWriter::Write(const BamRecordImpl& recordImpl)
{ d_->Write(internal::BamRecordMemory::GetRawData(recordImpl)); }
