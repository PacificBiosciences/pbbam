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
#include "pbbam/internal/BamRecordSort.h"
#include "pbbam/internal/MergeStrategy.h"
#include <algorithm>
#include <map>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace PacBio::BAM::staging;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

class ZmwQueryGroupIterator : public IBamFileGroupIterator
{
public:
    ZmwQueryGroupIterator(const std::vector<int>& zmwWhitelist,
                          const BamFile& file)
        : IBamFileGroupIterator(file)
    {
        std::vector<int> sortedZmws = zmwWhitelist;
        std::sort(sortedZmws.begin(), sortedZmws.end());
        for (int zmw : sortedZmws) {
            std::vector<int> zmwIndexList = { }; // PBI magic goes here for pbi.OffsetsForZmw(zmw);
            if (zmwIndexList.empty())
                continue;
            std::sort(zmwIndexList.begin(), zmwIndexList.end());
            zmwIndices_[zmw] = zmwIndexList;
        }
    }

public:
    bool GetNext(std::vector<BamRecord>& r) {

        if (zmwIndices_.empty())
            return false;

        const auto firstIter = zmwIndices_.cbegin();
        const std::vector<int>& indices = (*firstIter).second;
        for (int i : indices) {

            // seek & read
//            r.push_back( fileData_.records.at(i) );
        }
        zmwIndices_.erase(firstIter);
        return true;
    }

    bool InSameGroup(const BamRecord& lhs, const BamRecord& rhs) const {
        return lhs.HoleNumber() == rhs.HoleNumber();
    }

private:
    map<int, vector<int> > zmwIndices_;
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

ZmwGroupQuery::ZmwGroupQuery(const std::vector<int>& zmwWhitelist,
                             const DataSet& dataset)
    : IGroupQuery(dataset)
    , whitelist_(zmwWhitelist)
{
    mergeStrategy_.reset(new GroupMergeStrategy<ByZmw>(CreateIterators()));
}

ZmwGroupQuery::FileIterPtr ZmwGroupQuery::CreateIterator(const BamFile& file)
{ return FileIterPtr(new ZmwQueryGroupIterator(whitelist_, file)); }
