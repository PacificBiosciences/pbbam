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

#include "pbbam/PbiFilterZmwGroupQuery.h"
#include "CompositeBamReader.h"
#include <boost/optional.hpp>
#include <cassert>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace std;

struct PbiFilterZmwGroupQuery::PbiFilterZmwGroupQueryPrivate
{
    /// TODO: To ensure that BamRecords of a zmw are always stored 
    /// sequentially in bam files, such as in movie.subreads.bam.
public:
    PbiFilterZmwGroupQueryPrivate(const PbiFilter& filter, const DataSet& dataset)
        : reader_(new internal::PbiFilterCompositeBamReader<Compare::None>(filter, dataset))
        , nextRecord_(boost::none)
    { }

    bool GetNext(vector<BamRecord>& records)
    {
        records.clear();

        string movieName;
        uint32_t holeNumber;

        if (nextRecord_.is_initialized()) {
            BamRecord r = nextRecord_.get();
            movieName = r.MovieName();
            holeNumber = r.HoleNumber();
            records.push_back(std::move(r));
            nextRecord_ = boost::none;
        }

        BamRecord record;
        while (reader_->GetNext(record)) {
            if (records.empty()) {
                movieName = record.MovieName();
                holeNumber = record.HoleNumber();
                records.push_back(record);
            }
            else {
                assert(!records.empty());
                if (record.MovieName() == movieName and record.HoleNumber() == holeNumber)
                    records.push_back(record);
                else {
                    nextRecord_ = record;
                    return true;
                }
            }
        }
        return !records.empty();
    }

public:
    unique_ptr<internal::PbiFilterCompositeBamReader<Compare::None>> reader_;

    boost::optional<BamRecord> nextRecord_;
};

PbiFilterZmwGroupQuery::PbiFilterZmwGroupQuery(const DataSet& dataset)
    : internal::IGroupQuery()
    , d_(new PbiFilterZmwGroupQueryPrivate(  (dataset.Filters().Size() == 0) 
                                           ? (PbiFilter{PbiQueryLengthFilter{ 0 , Compare::GREATER_THAN } })
                                           : (PbiFilter::FromDataSet(dataset))
                                           , dataset))
{ }

PbiFilterZmwGroupQuery::PbiFilterZmwGroupQuery(const PbiFilter& filter, const DataSet& dataset)
    : internal::IGroupQuery()
    , d_(new PbiFilterZmwGroupQueryPrivate(filter, dataset))
{ }

PbiFilterZmwGroupQuery::~PbiFilterZmwGroupQuery(void) { }

bool PbiFilterZmwGroupQuery::GetNext(vector<BamRecord>& records)
{ return d_->GetNext(records); }
