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
/// \file BaiIndexedBamReader.cpp
/// \brief Implements the BaiIndexedBamReader class.
//
// Author: Derek Barnett

#include "pbbam/BaiIndexedBamReader.h"
#include "MemoryUtils.h"
using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

struct BaiIndexedBamReaderPrivate
{
public:
    BaiIndexedBamReaderPrivate(const BamFile& file,
                               const GenomicInterval& interval)
        : htsIndex_(nullptr)
        , htsIterator_(nullptr)
    {
        LoadIndex(file.Filename());
        Interval(file.Header(), interval);
    }

    void Interval(const BamHeader& header,
                  const GenomicInterval& interval)
    {
        htsIterator_.reset(nullptr);

        if (header.HasSequence(interval.Name())) {
            auto id = header.SequenceId(interval.Name());
            if (id >= 0 && static_cast<size_t>(id) < header.NumSequences()) {
                htsIterator_.reset(bam_itr_queryi(htsIndex_.get(),
                                                  id,
                                                  interval.Start(),
                                                  interval.Stop()));
            }
        }

        if (!htsIterator_)
            throw std::runtime_error("could not create iterator for requested region");
    }

    void LoadIndex(const string& fn)
    {
        htsIndex_.reset(bam_index_load(fn.c_str()));
        if (!htsIndex_)
            throw std::runtime_error("could not load BAI index data");
    }

    int ReadRawData(BGZF* bgzf, bam1_t* b)
    {
        assert(htsIterator_.get());
        return hts_itr_next(bgzf, htsIterator_.get(), b, nullptr);
    }

public:
    GenomicInterval interval_;
    std::unique_ptr<hts_idx_t, internal::HtslibIndexDeleter>    htsIndex_;
    std::unique_ptr<hts_itr_t, internal::HtslibIteratorDeleter> htsIterator_;
};

} // namespace internal
} // namespace BAM
} // namespace PacBio

BaiIndexedBamReader::BaiIndexedBamReader(const GenomicInterval& interval,
                                         const std::string& filename)
    : BaiIndexedBamReader(interval, BamFile(filename))
{ }

BaiIndexedBamReader::BaiIndexedBamReader(const GenomicInterval& interval,
                                         const BamFile& bamFile)
    : BamReader(bamFile)
    , d_(new BaiIndexedBamReaderPrivate(File(), interval))
{ }

BaiIndexedBamReader::BaiIndexedBamReader(const GenomicInterval& interval,
                                         BamFile&& bamFile)
    : BamReader(std::move(bamFile))
    , d_(new BaiIndexedBamReaderPrivate(File(), interval))
{ }

const GenomicInterval& BaiIndexedBamReader::Interval(void) const
{
    assert(d_);
    return d_->interval_;
}

int BaiIndexedBamReader::ReadRawData(BGZF* bgzf, bam1_t* b)
{
    assert(d_);
    return d_->ReadRawData(bgzf, b);
}

BaiIndexedBamReader& BaiIndexedBamReader::Interval(const GenomicInterval& interval)
{
    assert(d_);
    d_->Interval(Header(), interval);
    return *this;
}
