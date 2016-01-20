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
/// \file BamReader.cpp
/// \brief Implements the BamReader class.
//
// Author: Derek Barnett

#include "pbbam/BamReader.h"
#include "MemoryUtils.h"
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <cassert>
#include <cstdio>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

struct BamReaderPrivate
{
public:
    BamReaderPrivate(const BamFile& bamFile)
        : htsFile_(nullptr)
        , bamFile_(bamFile)
    {
        DoOpen();
    }

    BamReaderPrivate(BamFile&& bamFile)
        : htsFile_(nullptr)
        , bamFile_(std::move(bamFile))
    {
        DoOpen();
    }

    void DoOpen(void) {

        // fetch file pointer
        htsFile_.reset(sam_open(bamFile_.Filename().c_str(), "rb"));
        if (!htsFile_)
            throw std::runtime_error("could not open BAM file for reading");
    }

public:
    std::unique_ptr<samFile, internal::HtslibFileDeleter> htsFile_;
    BamFile bamFile_;
};

} // namespace internal
} // namespace BAM
} // namespace PacBio

BamReader::BamReader(const string& fn)
    : BamReader(BamFile(fn))
{ }

BamReader::BamReader(const BamFile& bamFile)
    : d_(new internal::BamReaderPrivate(bamFile))
{
    // skip header
    VirtualSeek(d_->bamFile_.FirstAlignmentOffset());
}

BamReader::BamReader(BamFile&& bamFile)
    : d_(new internal::BamReaderPrivate(std::move(bamFile)))
{
    // skip header
    VirtualSeek(d_->bamFile_.FirstAlignmentOffset());
}

BamReader::~BamReader(void) { }

BGZF* BamReader::Bgzf(void) const
{
    assert(d_);
    assert(d_->htsFile_);
    assert(d_->htsFile_->fp.bgzf);
    return d_->htsFile_->fp.bgzf;
}

const BamFile& BamReader::File(void) const
{
    assert(d_);
    return d_->bamFile_;
}

std::string BamReader::Filename(void) const
{
    assert(d_);
    return d_->bamFile_.Filename();
}

const BamHeader& BamReader::Header(void) const
{
    assert(d_);
    return d_->bamFile_.Header();
}

bool BamReader::GetNext(BamRecord& record)
{
    assert(Bgzf());
    assert(internal::BamRecordMemory::GetRawData(record).get());

    auto result = ReadRawData(Bgzf(), internal::BamRecordMemory::GetRawData(record).get());

    // success
    if (result >= 0) {
        internal::BamRecordMemory::UpdateRecordTags(record);
        record.header_ = Header();
        return true;
    }

    // EOF or end-of-data range (not an error)
    else if (result == -1)
        return false;

    // error corrupted file
    else {
        auto errorMsg = string{"corrupted BAM file: "};
        if (result == -2)
            errorMsg += "probably truncated";
        else if (result == -3)
            errorMsg += "could not read BAM record's' core data";
        else if (result == -4)
            errorMsg += "could not read BAM record's' variable-length data";
        else
            errorMsg += "unknown reason " + to_string(result);
        errorMsg += string{" ("};
        errorMsg += Filename();
        errorMsg += string{")"};
        throw std::runtime_error{errorMsg};
    }
}

int BamReader::ReadRawData(BGZF* bgzf, bam1_t* b)
{
    return bam_read1(bgzf, b);
}

void BamReader::VirtualSeek(int64_t virtualOffset)
{
    auto result = bgzf_seek(Bgzf(), virtualOffset, SEEK_SET);
    if (result != 0)
        throw std::runtime_error("Failed to seek in BAM file");
}

int64_t BamReader::VirtualTell(void) const
{
    return bgzf_tell(Bgzf());
}
