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
/// \file ZmwQuery.cpp
/// \brief Implements the ZmwQuery class.
//
// Author: Derek Barnett

#include "pbbam/ZmwGroupQuery.h"
#include "pbbam/BamRecord.h"
#include "pbbam/CompositeBamReader.h"
#include "pbbam/PbiFilterTypes.h"
#include "MemoryUtils.h"
#include <algorithm>
#include <deque>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace std;

struct ZmwGroupQuery::ZmwGroupQueryPrivate
{
    typedef PbiFilterCompositeBamReader<Compare::Zmw> ReaderType;
    typedef std::unique_ptr<ReaderType> ReaderPtr;

    ZmwGroupQueryPrivate(const std::vector<int32_t>& zmwWhitelist,
                         const DataSet& dataset)
        : whitelist_(zmwWhitelist.cbegin(), zmwWhitelist.cend())
        , reader_(nullptr)
    {
        std::sort(whitelist_.begin(), whitelist_.end());
        whitelist_.erase(std::unique(whitelist_.begin(),
                                     whitelist_.end()),
                         whitelist_.end());

        if (!whitelist_.empty()) {
            reader_ = ReaderPtr(new ReaderType(PbiZmwFilter{whitelist_.front()}, dataset));
            whitelist_.pop_front();
        }
    }

    bool GetNext(std::vector<BamRecord>& records)
    {
        records.clear();
        if (!reader_)
            return false;

        // get all records matching ZMW
        BamRecord r;
        while (reader_->GetNext(r))
            records.push_back(r);

        // set next ZMW (if any left)
        if (!whitelist_.empty()) {
            reader_->Filter(PbiZmwFilter{whitelist_.front()});
            whitelist_.pop_front();
        }

        // otherwise destroy reader, next iteration will return false
        else
            reader_.reset(nullptr);

        return true;
    }

    std::deque<int32_t> whitelist_;
    ReaderPtr reader_;
};

ZmwGroupQuery::ZmwGroupQuery(const std::vector<int32_t>& zmwWhitelist,
                             const DataSet& dataset)
    : internal::IGroupQuery()
    , d_(new ZmwGroupQueryPrivate(zmwWhitelist, dataset))
{ }

ZmwGroupQuery::~ZmwGroupQuery(void) { }

bool ZmwGroupQuery::GetNext(std::vector<BamRecord>& records)
{ return d_->GetNext(records); }
