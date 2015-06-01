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
#include "pbbam/internal/BamRecordSort.h"
#include "pbbam/internal/MergeStrategy.h"
#include <algorithm>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace PacBio::BAM::staging;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

class ZmwQueryIterator : public IBamFileIterator
{
public:
    ZmwQueryIterator(const std::vector<int>& zmwWhitelist,
                     const BamFile& bamFile)
        : internal::IBamFileIterator(bamFile)
        , currentWhitelistIndex_(0)
    {
        std::vector<int> sortedZmws = zmwWhitelist;
        std::sort(sortedZmws.begin(), sortedZmws.end());
        for (int zmw : sortedZmws) {
            (void)zmw;
            std::vector<int> zmwIndices = { }; // PBI magic goes here for pbi.OffsetsForZmw(zmw);
            std::sort(zmwIndices.begin(), zmwIndices.end());
            for (int index : zmwIndices)
                whitelistRecordIndices_.push_back(index);
        }
    }

public:
    bool GetNext(BamRecord& r){
        if (currentWhitelistIndex_ >= whitelistRecordIndices_.size())
            return false;

        // this is where we seek & read
        // r = fileData_.records.at(whitelistRecordIndices_.at(currentWhitelistIndex_));

        ++currentWhitelistIndex_;
        return true;
    }

private:
    std::vector<int> whitelistRecordIndices_;
    size_t currentWhitelistIndex_;
};

} // namespace internal
} // namespace BAM
} // namespace PacBio

ZmwQuery::ZmwQuery(const std::vector<int>& zmwWhitelist,
                   const DataSet& dataset)
    : internal::IQuery(dataset)
    , whitelist_(zmwWhitelist)
{
    // not yet fully implemented
    throw std::exception();

//    mergeStrategy_.reset(new MergeStrategy<ByZmw>(CreateIterators()));
}

ZmwQuery::FileIterPtr ZmwQuery::CreateIterator(const BamFile& bamFile)
{ return FileIterPtr(new ZmwQueryIterator(whitelist_, bamFile)); }
