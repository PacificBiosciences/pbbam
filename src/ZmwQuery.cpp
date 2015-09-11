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

#include "pbbam/ZmwQuery.h"
#include "pbbam/PbiIndex.h"
#include "pbbam/internal/BamRecordSort.h"
#include "pbbam/internal/MergeStrategy.h"
#include "MemoryUtils.h"
#include <htslib/bgzf.h>
#include <htslib/sam.h>
#include <algorithm>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
//using namespace PacBio::BAM::staging;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

class ZmwQueryIterator : public IBamFileIterator
{
public:
    ZmwQueryIterator(const std::vector<int32_t>& zmwWhitelist,
                     const BamFile& bamFile)
        : internal::IBamFileIterator(bamFile)
        , currentBlockReadCount_(0)
        , htsFile_(nullptr)
        , htsHeader_(nullptr)
    {
        // init BAM file for reading
        htsFile_.reset(sam_open(bamFile.Filename().c_str(), "rb"));
        if (!htsFile_)
            throw std::runtime_error("could not open BAM file for reading");

        htsHeader_.reset(sam_hdr_read(htsFile_.get()));
        if (!htsHeader_)
            throw std::runtime_error("could not read BAM header data");

        // open index & query for ZMWs
        PbiIndex index(bamFile.PacBioIndexFilename());
        blocks_ = index.Lookup(ZmwIndexMultiRequest(zmwWhitelist));
    }

public:
    bool GetNext(BamRecord& r){

        // no data to fetch, return false
        if (blocks_.empty())
            return false;

        // maybe seek to block
        if (currentBlockReadCount_ == 0) {
            const int seekResult = bgzf_seek(htsFile_.get()->fp.bgzf, blocks_.at(0).virtualOffset_, SEEK_SET);
            if (seekResult == -1)
                throw std::runtime_error("could not seek in BAM file");
        }

        // read next record
//        r = BamRecord(fileData_.Header());
        const int readResult = sam_read1(htsFile_.get(),
                                         htsHeader_.get(),
                                         internal::BamRecordMemory::GetRawData(r).get());
        internal::BamRecordMemory::UpdateRecordTags(r);
//        r.header_ = fileData_.Header();
        r.header_ = header_;

        // update counters
        ++currentBlockReadCount_;
        if (currentBlockReadCount_ == blocks_.at(0).numReads_) {
            blocks_.pop_front();
            currentBlockReadCount_ = 0;
        }

        // return result of reading
        if (readResult >= 0)           // success
            return true;
        else if (readResult == -1)     // normal EOF
            return false;
        else                           // error (truncated file, etc)
            throw std::runtime_error("corrupted file, may be truncated");
    }

private:
    IndexResultBlocks blocks_;
    size_t currentBlockReadCount_;
    unique_ptr<samFile,   internal::HtslibFileDeleter>   htsFile_;
    unique_ptr<bam_hdr_t, internal::HtslibHeaderDeleter> htsHeader_;
};

} // namespace internal
} // namespace BAM
} // namespace PacBio

ZmwQuery::ZmwQuery(const std::vector<int32_t> &zmwWhitelist,
                   const DataSet& dataset)
    : internal::IQuery(dataset)
    , whitelist_(zmwWhitelist)
{
    mergeStrategy_.reset(new MergeStrategy<ByZmw>(CreateIterators()));
}

ZmwQuery::FileIterPtr ZmwQuery::CreateIterator(const BamFile& bamFile)
{ return FileIterPtr(new ZmwQueryIterator(whitelist_, bamFile)); }
