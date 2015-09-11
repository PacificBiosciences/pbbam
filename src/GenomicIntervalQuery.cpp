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

#include "pbbam/GenomicIntervalQuery.h"
#include "pbbam/BamFile.h"
#include "pbbam/internal/BamRecordSort.h"
#include "pbbam/internal/MergeStrategy.h"
#include "AssertUtils.h"
#include "MemoryUtils.h"
#include <cassert>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

class GenomicIntervalIterator : public internal::IBamFileIterator
{
public:
    GenomicIntervalIterator(const GenomicInterval& interval,
                            const BamFile& bamFile)
        : internal::IBamFileIterator(bamFile)
        , interval_(interval)
    {
        // open file
        htsFile_.reset(sam_open(bamFile.Filename().c_str(), "rb"));
        if (!htsFile_)
            throw std::runtime_error("could not open BAM file for reading");

        // load header info
        htsHeader_.reset(sam_hdr_read(htsFile_.get()));
        if (!htsHeader_)
            throw std::runtime_error("could not read BAM header data");

        // open index
        htsIndex_.reset(bam_index_load(bamFile.Filename().c_str()));
        if (!htsIndex_)
            throw std::runtime_error("could not load BAI index data");

        // initialize iterator
        if (bamFile.Header().HasSequence(interval_.Name())) {
            const int id = bamFile.ReferenceId(interval_.Name());
            if (id >= 0 && id < htsHeader_->n_targets) {
                htsIterator_.reset(sam_itr_queryi(htsIndex_.get(),
                                                  id,
                                                  interval.Start(),
                                                  interval.Stop()));
            }
        }
        if (!htsIterator_)
            throw std::runtime_error("could not create iterator for requested region");
    }

public:
    bool GetNext(BamRecord& record) {

        assert(htsFile_);
        assert(htsIterator_);
        const int result = sam_itr_next(htsFile_.get(),
                                        htsIterator_.get(),
                                        internal::BamRecordMemory::GetRawData(record).get());
        internal::BamRecordMemory::UpdateRecordTags(record);
        record.header_ = header_;

        // success
        if (result >= 0)
            return true;

        // normal EOF
        else if (result == -1)
            return false;

        // error (truncated file, etc)
        else
            throw std::runtime_error("corrupted file, may be truncated");
    }

private:
    GenomicInterval interval_;
    unique_ptr<samFile,   internal::HtslibFileDeleter>     htsFile_;
    unique_ptr<bam_hdr_t, internal::HtslibHeaderDeleter>   htsHeader_;
    unique_ptr<hts_idx_t, internal::HtslibIndexDeleter>    htsIndex_;
    unique_ptr<hts_itr_t, internal::HtslibIteratorDeleter> htsIterator_;
};

GenomicIntervalQuery::GenomicIntervalQuery(const GenomicInterval& interval,
                                           const DataSet& dataset)
    : internal::IQuery(dataset)
    , interval_(interval)
{
    Interval(interval_);
}

GenomicIntervalQuery::FileIterPtr GenomicIntervalQuery::CreateIterator(const BamFile& bamFile)
{ return FileIterPtr(new GenomicIntervalIterator(interval_, bamFile)); }

GenomicIntervalQuery& GenomicIntervalQuery::Interval(const GenomicInterval& interval)
{
    interval_ = interval;
    // check files
    // if SO all coordinate
    // else if SO all queryname
    // else SO unsorted/unknown
    mergeStrategy_.reset(new internal::MergeStrategy<internal::ByPosition>(CreateIterators()));
    return *this;
}

GenomicInterval GenomicIntervalQuery::Interval(void) const
{ return interval_; }
