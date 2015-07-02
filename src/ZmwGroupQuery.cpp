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

#include "pbbam/ZmwGroupQuery.h"
#include "pbbam/PbiIndex.h"
#include "pbbam/internal/BamRecordSort.h"
#include "pbbam/internal/MergeStrategy.h"
#include "MemoryUtils.h"
#include <algorithm>
#include <map>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
//using namespace PacBio::BAM::staging;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

class ZmwQueryGroupIterator : public IBamFileGroupIterator
{
public:
    ZmwQueryGroupIterator(const std::vector<int32_t>& zmwWhitelist,
                          const BamFile& file)
        : IBamFileGroupIterator(file)
    {
        // init BAM file for reading
        htsFile_.reset(sam_open(file.Filename().c_str(), "rb"));
        if (!htsFile_)
            throw std::runtime_error("could not open BAM file for reading");

        htsHeader_.reset(sam_hdr_read(htsFile_.get()));
        if (!htsHeader_)
            throw std::runtime_error("could not read BAM header data");

        // open index & query for ZMWs
        PbiIndex index(file.PacBioIndexFilename());
        for (int32_t zmw : zmwWhitelist)
            zmwGroups_[zmw] = index.Lookup(ZmwIndexRequest(zmw));
    }

public:
    bool GetNext(std::vector<BamRecord>& r)
    {
        r.clear();
        if (zmwGroups_.empty())
            return false;

        BamRecord record(header_);
        const IndexResultBlocks& blocks = zmwGroups_.cbegin()->second;
        for (const IndexResultBlock& block : blocks) {

            // seek to first record in block
            const int seekResult = bgzf_seek(htsFile_.get()->fp.bgzf, block.virtualOffset_, SEEK_SET);
            if (seekResult == -1)
                throw std::runtime_error("could not seek in BAM file");

            // read block records
            for (size_t i = 0; i < block.numReads_; ++i) {
                const int readResult = sam_read1(htsFile_.get(),
                                                 htsHeader_.get(),
                                                 internal::BamRecordMemory::GetRawData(record).get());
//                record.header_ = fileData_.Header();

                if (readResult >= 0)           // success
                    r.push_back(record);
                else if (readResult == -1)     // normal EOF
                    break;
                else                           // error (truncated file, etc)
                    throw std::runtime_error("corrupted file, may be truncated");
            }
        }

        // pop zmw info & return success
        zmwGroups_.erase(zmwGroups_.begin());
        return !r.empty();
    }

    bool InSameGroup(const BamRecord& lhs, const BamRecord& rhs) const
    { return lhs.HoleNumber() == rhs.HoleNumber(); }

private:
    unique_ptr<samFile,   internal::HtslibFileDeleter>   htsFile_;
    unique_ptr<bam_hdr_t, internal::HtslibHeaderDeleter> htsHeader_;
    map<int32_t, IndexResultBlocks> zmwGroups_;
};

} // namespace internal
} // namespace BAM
} // namespace PacBio

ZmwGroupQuery::ZmwGroupQuery(const DataSet& dataset)
    : IGroupQuery(dataset)
    , whitelist_(/* all dataset ZMWs */)
{
    mergeStrategy_.reset(new GroupMergeStrategy<ByZmw>(CreateIterators()));
}

ZmwGroupQuery::ZmwGroupQuery(const std::vector<int32_t>& zmwWhitelist,
                             const DataSet& dataset)
    : IGroupQuery(dataset)
    , whitelist_(zmwWhitelist)
{
    mergeStrategy_.reset(new GroupMergeStrategy<ByZmw>(CreateIterators()));
}

ZmwGroupQuery::FileIterPtr ZmwGroupQuery::CreateIterator(const BamFile& file)
{ return FileIterPtr(new ZmwQueryGroupIterator(whitelist_, file)); }
