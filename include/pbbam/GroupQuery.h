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

// Author: Yuan Li

#ifndef _GROUP_QUERY_H_
#define _GROUP_QUERY_H_
#include "GroupQueryBase.h"
#include <htslib/sam.h>
#include <vector>

namespace PacBio {
namespace BAM {

class PBBAM_EXPORT SequentialGroupQueryBase: public GroupQueryBase
{
public:
    SequentialGroupQueryBase(const BamFile & bamFile);

protected:
    virtual bool InSameGroup(const BamRecord & record, const BamRecord & another) = 0;
    bool GetNext(std::vector<BamRecord> & records);
    PBBAM_SHARED_PTR<samFile>   htsFile_;
    PBBAM_SHARED_PTR<bam_hdr_t> htsHeader_;
    BamRecord nextRecord_;
};

class PBBAM_EXPORT ZmwQuery: public SequentialGroupQueryBase
{
public: 
    ZmwQuery(const BamFile & bamFile)
    : SequentialGroupQueryBase(bamFile) { }

private:
    bool InSameGroup(const BamRecord & record, const BamRecord & another) {
        return (record.MovieName() == another.MovieName() and 
                record.HoleNumber() == another.HoleNumber());
    }
};

class PBBAM_EXPORT QNameQuery: public SequentialGroupQueryBase
{
public:
    QNameQuery(const BamFile & bamFile) 
    : SequentialGroupQueryBase(bamFile) { }

private:
    bool InSameGroup(const BamRecord & record, const BamRecord & another) {
        return (record.Impl().Name() == another.Impl().Name());
    }
};

} // namespace BAM
} // namespace PacBio

#endif // _SEQUENTIAL_GROUP_QUERY_BASE_H_
