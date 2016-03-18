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
#include "pbbam/Validator.h"
#include "AssertUtils.h"
#include "FileProducer.h"
#include "MemoryUtils.h"
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <iostream>
#include <thread>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

class BamWriterPrivate : public internal::FileProducer
{
public:
    BamWriterPrivate(const std::string& filename,
                     const PBBAM_SHARED_PTR<bam_hdr_t> rawHeader,
                     const BamWriter::CompressionLevel compressionLevel,
                     const size_t numThreads,
                     const BamWriter::BinCalculationMode binCalculationMode);

public:
    void Write(const BamRecord& record);
    void Write(const BamRecord& record, int64_t* vOffset);
    void Write(const BamRecordImpl& recordImpl);

public:
    bool calculateBins_;
    std::unique_ptr<samFile, internal::HtslibFileDeleter> file_;
    PBBAM_SHARED_PTR<bam_hdr_t> header_;
};

BamWriterPrivate::BamWriterPrivate(const string& filename,
                                   const PBBAM_SHARED_PTR<bam_hdr_t> rawHeader,
                                   const BamWriter::CompressionLevel compressionLevel,
                                   const size_t numThreads,
                                   const BamWriter::BinCalculationMode binCalculationMode)
    : internal::FileProducer(filename)
    , calculateBins_(binCalculationMode == BamWriter::BinCalculation_ON)
    , file_(nullptr)
    , header_(rawHeader)
{
    if (!header_)
        throw std::runtime_error("null header");

    // open file
    const string& usingFilename = TempFilename();
    const string& mode = string("wb") + to_string(static_cast<int>(compressionLevel));
    file_.reset(sam_open(usingFilename.c_str(), mode.c_str()));
    if (!file_)
        throw std::runtime_error("could not open file for writing");

    // if no explicit thread count given, attempt built-in check
    size_t actualNumThreads = numThreads;
    if (actualNumThreads == 0) {
        actualNumThreads = thread::hardware_concurrency();

        // if still unknown, default to single-threaded
        if (actualNumThreads == 0)
            actualNumThreads = 1;
    }

    // if multithreading requested, enable it
    if (actualNumThreads > 1)
        hts_set_threads(file_.get(), actualNumThreads);

    // write header
    const int ret = sam_hdr_write(file_.get(), header_.get());
    if (ret != 0)
        throw std::runtime_error("could not write header");
}

void BamWriterPrivate::Write(const BamRecord& record)
{
#if PBBAM_AUTOVALIDATE
    Validator::Validate(record);
#endif

    const auto rawRecord = internal::BamRecordMemory::GetRawData(record);

    // (probably) store bins
    // min_shift=14 & n_lvls=5 are BAM "magic numbers"
    if (calculateBins_)
        rawRecord->core.bin = hts_reg2bin(rawRecord->core.pos, bam_endpos(rawRecord.get()), 14, 5);

    // write record to file
    const int ret = sam_write1(file_.get(), header_.get(), rawRecord.get());
    if (ret <= 0)
        throw std::runtime_error("could not write record");
}

void BamWriterPrivate::Write(const BamRecord& record, int64_t* vOffset)
{
    BGZF* bgzf = file_.get()->fp.bgzf;
    assert(bgzf);
    assert(vOffset);

    // ensure offsets up-to-date
    bgzf_flush(bgzf);

    // capture virtual offset where we’re about to write
    const off_t rawTell = htell(bgzf->fp);
    const int length = bgzf->block_offset;
    *vOffset = (rawTell << 16) | length ;

    // now write data
    Write(record);
}

inline void BamWriterPrivate::Write(const BamRecordImpl& recordImpl)
{ Write(BamRecord(recordImpl)); }

} // namespace internal
} // namespace BAM
} // namespace PacBio

BamWriter::BamWriter(const std::string& filename,
                     const BamHeader& header,
                     const BamWriter::CompressionLevel compressionLevel,
                     const size_t numThreads,
                     const BinCalculationMode binCalculationMode)
    : d_(nullptr)
{
#if PBBAM_AUTOVALIDATE
    Validator::Validate(header);
#endif
    d_.reset(new internal::BamWriterPrivate{ filename,
                                             internal::BamHeaderMemory::MakeRawHeader(header),
                                             compressionLevel,
                                             numThreads,
                                             binCalculationMode
                                           });
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
        throw std::runtime_error("could not flush output buffer contents");
}

void BamWriter::Write(const BamRecord& record)
{ d_->Write(record); }

void BamWriter::Write(const BamRecord& record, int64_t* vOffset)
{ d_->Write(record, vOffset); }

void BamWriter::Write(const BamRecordImpl& recordImpl)
{ d_->Write(recordImpl); }
