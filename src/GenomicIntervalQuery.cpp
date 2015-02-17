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

#include "pbbam/GenomicIntervalQuery.h"
#include "pbbam/BamFile.h"
#include "AssertUtils.h"
#include "MemoryUtils.h"
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

GenomicIntervalQuery::GenomicIntervalQuery(const GenomicInterval& interval,
                                           const BamFile& file)
    : QueryBase()
{
    if (InitFile(file))
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
    if (error_ == GenomicIntervalQuery::NoError && iterator_) {
        const int result = sam_itr_next(file_.get(), iterator_.get(), record.RawData().get());
        if (result >= 0)
            return true;
        else if (result < -1) {
            // TODO: determine & report error (truncated BAM file or malformed index, etc.)
        }
    }
    return false;
}

bool GenomicIntervalQuery::InitFile(const BamFile& file) {

    error_ = GenomicIntervalQuery::NoError;

    // open file
    file_.reset(sam_open(file.Filename().c_str(), "rb"), internal::HtslibFileDeleter());
    if (!file_) {
        error_ = GenomicIntervalQuery::FileOpenError;
        return false;
    }

    header_.reset(sam_hdr_read(file_.get()), internal::HtslibHeaderDeleter());
    if (!header_) {
        error_ = GenomicIntervalQuery::FileMetadataError;
        return false;
    }

    // open index
    index_.reset(bam_index_load(file.Filename().c_str()), internal::HtslibIndexDeleter());
    if (!index_) {
        error_ = GenomicIntervalQuery::IndexFileOpenError;
        return false;
    }

    // if we get here, return success
    return true;
}

GenomicInterval GenomicIntervalQuery::Interval(void) const
{
    return interval_;
}

GenomicIntervalQuery& GenomicIntervalQuery::Interval(const GenomicInterval& interval)
{
    // if file-related error, or missing data - then setting a new interval
    // can't help anything. just get out of here
    if (!file_ || !header_ || !index_ ||
        error_ == GenomicIntervalQuery::FileOpenError ||
        error_ == GenomicIntervalQuery::FileMetadataError ||
        error_ == GenomicIntervalQuery::IndexFileOpenError ||
        error_ == GenomicIntervalQuery::IndexFileMetadataError)
    {
        iterator_.reset();
        return *this;
    }

    // otherwise, attempt to get an iterator for this interval
    error_ = GenomicIntervalQuery::InitializeQueryError;
    if (interval.Id() >= 0 && interval.Id() < header_->n_targets) {
        iterator_.reset(sam_itr_queryi(index_.get(),
                                       interval.Id(),
                                       interval.Start(),
                                       interval.Stop()),
                        internal::HtslibIteratorDeleter());
        // if successful, clear error
        if(iterator_)
            error_ = GenomicIntervalQuery::NoError;
    }

    return *this;
}
