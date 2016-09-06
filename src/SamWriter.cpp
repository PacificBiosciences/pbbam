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

#include "pbbam/SamWriter.h"
#include "pbbam/Validator.h"
#include "FileProducer.h"
#include "MemoryUtils.h"
#include <htslib/hfile.h>
#include <htslib/sam.h>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

class SamWriterPrivate : public internal::FileProducer
{
public:
    SamWriterPrivate(const std::string& filename,
                      const PBBAM_SHARED_PTR<bam_hdr_t> rawHeader)
        : internal::FileProducer(filename)
        , file_(nullptr)
        , header_(rawHeader)
    {
        if (!header_)
            throw std::runtime_error("null header");

        // open file
        const string& usingFilename = TempFilename();
        const string& mode = string("w");
        file_.reset(sam_open(usingFilename.c_str(), mode.c_str()));
        if (!file_)
            throw std::runtime_error("could not open file for writing");

        // write header
        const int ret = sam_hdr_write(file_.get(), header_.get());
        if (ret != 0)
            throw std::runtime_error("could not write header");
    }

    void TryFlush(void);
    void Write(const BamRecord& record);

private:
    std::unique_ptr<samFile, internal::HtslibFileDeleter> file_;
    PBBAM_SHARED_PTR<bam_hdr_t> header_;
};

void SamWriterPrivate::TryFlush(void)
{
    const auto ret = file_.get()->fp.hfile;
    if (ret != 0)
        throw std::runtime_error("could not flush output buffer contents");
}

void SamWriterPrivate::Write(const BamRecord& record)
{
#if PBBAM_AUTOVALIDATE
    Validator::Validate(record);
#endif

    const auto rawRecord = internal::BamRecordMemory::GetRawData(record);

    // store bin number
    // min_shift=14 & n_lvls=5 are SAM/BAM "magic numbers"
    rawRecord->core.bin = hts_reg2bin(rawRecord->core.pos,
                                      bam_endpos(rawRecord.get()), 14, 5);

    // write record to file
    const int ret = sam_write1(file_.get(), header_.get(), rawRecord.get());
    if (ret <= 0)
        throw std::runtime_error("could not write record");
}

} // namespace internal
} // namespace BAM
} // namespace PacBio

SamWriter::SamWriter(const string& filename, const BamHeader& header)
    : IRecordWriter()
    , d_(nullptr)
{
#if PBBAM_AUTOVALIDATE
    Validator::Validate(header);
#endif
    d_.reset(new internal::SamWriterPrivate{ filename,
                                             internal::BamHeaderMemory::MakeRawHeader(header)
                                           });
}

SamWriter::~SamWriter(void) { }

void SamWriter::TryFlush(void)
{
    d_->TryFlush();
}

void SamWriter::Write(const BamRecord& record)
{
    d_->Write(record);
}

void SamWriter::Write(const BamRecordImpl& recordImpl)
{
    d_->Write( BamRecord{recordImpl} );
}
