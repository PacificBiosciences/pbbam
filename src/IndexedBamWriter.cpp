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
/// \file IndexedBamWriter.cpp
/// \brief Implements the IndexedBamWriter class
//
// Author: Derek Barnett

#include "pbbam/IndexedBamWriter.h"

#include <cassert>
#include <cstdint>
#include <stdexcept>

#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>

#include "pbbam/BamHeader.h"
#include "pbbam/BamRecord.h"
#include "pbbam/PbiBuilder.h"
#include "pbbam/Unused.h"
#include "pbbam/Validator.h"

#include "FileProducer.h"
#include "MemoryUtils.h"

namespace PacBio {
namespace BAM {
namespace internal {

class IndexedBamWriterPrivate : public internal::FileProducer
{
public:
    IndexedBamWriterPrivate(const std::string& outputFilename,
                            const std::shared_ptr<bam_hdr_t> rawHeader)
        : internal::FileProducer{outputFilename}
        , file_{nullptr}
        , header_{rawHeader}
        , builder_{outputFilename + ".pbi"}
        , previousBlockAddress_{0}
    {
        if (!header_) throw std::runtime_error("null header");

        // open file
        const auto& usingFilename = TempFilename();
        file_.reset(sam_open(usingFilename.c_str(), "wb"));
        if (!file_) throw std::runtime_error("could not open file for writing");

        // write header
        const auto ret = sam_hdr_write(file_.get(), header_.get());
        if (ret != 0) throw std::runtime_error("could not write header");

        // store first alignment block
        previousBlockAddress_ = file_.get()->fp.bgzf->block_address;
    }

    ~IndexedBamWriterPrivate()
    {
        // ensure last remaining bits are flushed to file
        const auto ret = bgzf_flush(file_.get()->fp.bgzf);
        UNUSED(ret);
    }

public:
    void Write(const BamRecord& record)
    {
#if PBBAM_AUTOVALIDATE
        Validator::Validate(record);
#endif
        const auto rawRecord = internal::BamRecordMemory::GetRawData(record);
        BGZF* bgzf = file_.get()->fp.bgzf;
        assert(bgzf);

        //
        // Fetch record's start offset.
        //
        // If we're still in the same block from the last record written, we
        // need to flush to get the proper offset. This will
        //
        if (bgzf->block_address == previousBlockAddress_) bgzf_flush(bgzf);
        const int64_t vOffset = bgzf_tell(bgzf);

        // update bin
        rawRecord->core.bin = hts_reg2bin(rawRecord->core.pos, bam_endpos(rawRecord.get()), 14, 5);

        // write record to file & PBI builder
        const auto ret = sam_write1(file_.get(), header_.get(), rawRecord.get());
        if (ret <= 0) throw std::runtime_error("could not write record");
        builder_.AddRecord(record, vOffset);

        // update block address
        previousBlockAddress_ = bgzf->block_address;
    }

public:
    std::unique_ptr<samFile, internal::HtslibFileDeleter> file_;
    std::shared_ptr<bam_hdr_t> header_;
    PbiBuilder builder_;
    int64_t previousBlockAddress_;
};

}  // namespace internal

IndexedBamWriter::IndexedBamWriter(const std::string& outputFilename, const BamHeader& header)
    : IRecordWriter{}, d_{nullptr}
{
#if PBBAM_AUTOVALIDATE
    Validator::Validate(header);
#endif
    d_ = std::make_unique<internal::IndexedBamWriterPrivate>(
        outputFilename, internal::BamHeaderMemory::MakeRawHeader(header));
}

IndexedBamWriter::~IndexedBamWriter() {}

void IndexedBamWriter::TryFlush()
{
    const auto ret = bgzf_flush(d_->file_.get()->fp.bgzf);
    if (ret != 0) throw std::runtime_error("could not flush output buffer contents");
}

void IndexedBamWriter::Write(const BamRecord& record) { d_->Write(record); }

void IndexedBamWriter::Write(const BamRecordImpl& record) { d_->Write(BamRecord{record}); }

}  // namespace BAM
}  // namespace PacBio
