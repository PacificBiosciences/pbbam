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

GenomicIntervalQuery::GenomicIntervalQuery(const GenomicInterval& interval,
                                           const BamFile& file)
    : QueryBase(file)
    , htsIterator_(nullptr)
{
    // open file
    htsFile_.reset(sam_open(file.Filename().c_str(), "rb"), internal::HtslibFileDeleter());
    if (!htsFile_)
        throw std::exception();

    htsHeader_.reset(sam_hdr_read(htsFile_.get()), internal::HtslibHeaderDeleter());
    if (!htsHeader_)
        throw std::exception();

    // open index
    htsIndex_.reset(bam_index_load(file.Filename().c_str()), internal::HtslibIndexDeleter());
    if (!htsIndex_)
        throw std::exception();

    // initialize interval
    Interval(interval);
}

//GenomicIntervalQuery::GenomicIntervalQuery(const string& zeroBasedRegion,
//                                           const BamFile& file)
//    : QueryBase()
//{
//    if (InitFile(file))
//        Interval(zeroBasedRegion);
//}

bool GenomicIntervalQuery::GetNext(BamRecord& record)
{
    if (htsIterator_) {
        const int result = sam_itr_next(htsFile_.get(),
                                        htsIterator_.get(),
                                        internal::BamRecordMemory::GetRawData(record).get());

        // success
        if (result >= 0)
            return true;

        // normal EOF
        else if (result == -1)
            return false;

        // error (truncated file, etc)
        else
            throw std::exception();
    }

    // no iterator set
    return false;
}

GenomicInterval GenomicIntervalQuery::Interval(void) const
{ return interval_; }

GenomicIntervalQuery& GenomicIntervalQuery::Interval(const GenomicInterval& interval)
{
    // if file-related error, or missing data - setting a new interval
    // can't help anything. just get out of here
    if (!htsFile_ || !htsHeader_ || !htsIndex_) {
        throw std::exception();
    }

    // ensure clean slate
    htsIterator_.reset();

    // lookup ID for reference name
    if (file_.Header().HasSequence(interval.Name())) {
        const int id = file_.ReferenceId(interval.Name());
        if (id >= 0 && id < htsHeader_->n_targets) {

            // get iterator for interval
            htsIterator_.reset(sam_itr_queryi(htsIndex_.get(),
                                              id,
                                              interval.Start(),
                                              interval.Stop()),
                               internal::HtslibIteratorDeleter());
        }
    }

    return *this;
}

namespace PacBio {
namespace BAM {
namespace staging {

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
            throw std::exception();

        // load header info
        htsHeader_.reset(sam_hdr_read(htsFile_.get()));
        if (!htsHeader_)
            throw std::exception();

        // open index
        htsIndex_.reset(bam_index_load(bamFile.Filename().c_str()));
        if (!htsIndex_)
            throw std::exception();

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
            throw std::exception();
    }

public:
    bool GetNext(BamRecord& record) {

        assert(htsFile_);
        assert(htsIterator_);
        const int result = sam_itr_next(htsFile_.get(),
                                        htsIterator_.get(),
                                        internal::BamRecordMemory::GetRawData(record).get());
        // success
        if (result >= 0)
            return true;

        // normal EOF
        else if (result == -1)
            return false;

        // error (truncated file, etc)
        else
            throw std::exception();
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

} // namespace staging
} // namespace BAM
} // namespace PacBio
