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

#ifndef PBIFILTER_ZMWGROUPQUERY_H 
#define PBIFILTER_ZMWGROUPQUERY_H

#include "pbbam/internal/QueryBase.h"
#include "pbbam/PbiFilter.h"
#include <memory>

namespace PacBio {
namespace BAM {

/// This class operates on a name-sorted BAM file, with each iteration of the query
/// returning each contiguous block of records that share a name.
///
/// \note Iterate over zmws, return vector of subreads of a zmw each time.
///
class PBBAM_EXPORT PbiFilterZmwGroupQuery : public internal::IGroupQuery
{
public:
    PbiFilterZmwGroupQuery(const DataSet& dataset);
    PbiFilterZmwGroupQuery(const PbiFilter& filter, const DataSet& dataset);
    ~PbiFilterZmwGroupQuery(void);

public:
    bool GetNext(std::vector<BamRecord>& records);

private:
    struct PbiFilterZmwGroupQueryPrivate;
    std::unique_ptr<PbiFilterZmwGroupQueryPrivate> d_;
};

} // namespace BAM
} // namespace PacBio

#endif // PBIFILTER_ZMWGROUPQUERY_H
