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

#ifndef COMPOSITEBAMREADER_H
#define COMPOSITEBAMREADER_H

#include "pbbam/BamFile.h"
#include "pbbam/BamHeader.h"
#include "pbbam/BamRecord.h"
#include "pbbam/DataSet.h"
#include "pbbam/GenomicInterval.h"
#include "BaiIndexedBamReader.h"
#include "BamReader.h"
#include "PbiIndexedBamReader.h"

#include <boost/optional.hpp>

#include <deque>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {
namespace internal {

struct CompositeMergeItem
{
public:
    std::unique_ptr<BamReader> reader;
    BamRecord record;

public:
    CompositeMergeItem(std::unique_ptr<BamReader>&& rdr);
    CompositeMergeItem(std::unique_ptr<BamReader>&& rdr, BamRecord&& rec);
    CompositeMergeItem(CompositeMergeItem&& other);
    CompositeMergeItem& operator=(CompositeMergeItem&& other);
    ~CompositeMergeItem(void);
};

template<typename CompareType>
struct CompositeMergeItemSorter : public std::function<bool(const CompositeMergeItem&,
                                                            const CompositeMergeItem&)>
{
    bool operator()(const CompositeMergeItem& lhs,
                    const CompositeMergeItem& rhs);
};

class GenomicIntervalCompositeBamReader
{
public:
    GenomicIntervalCompositeBamReader(const GenomicInterval& interval, const std::vector<BamFile>& bamFiles);
    GenomicIntervalCompositeBamReader(const GenomicInterval& interval, std::vector<BamFile>&& bamFiles);
    GenomicIntervalCompositeBamReader(const GenomicInterval& interval, const DataSet& dataset);

public:
    bool GetNext(BamRecord& record);

public:
    GenomicIntervalCompositeBamReader& Interval(const GenomicInterval& interval);
    const GenomicInterval& Interval(void) const;

private:
    void UpdateSort(void);

private:
    GenomicInterval interval_;
    std::deque<CompositeMergeItem> mergeItems_;
    std::vector<std::string> filenames_;
};

template<typename OrderByType>
class PbiFilterCompositeBamReader
{
public:
    typedef CompositeMergeItem                              value_type;
    typedef CompositeMergeItemSorter<OrderByType>           merge_sorter_type;
    typedef std::deque<value_type>                  container_type;
    typedef typename container_type::iterator       iterator;
    typedef typename container_type::const_iterator const_iterator;

public:
    PbiFilterCompositeBamReader(const PbiFilter& filter, const std::vector<BamFile>& bamFiles);
    PbiFilterCompositeBamReader(const PbiFilter& filter, std::vector<BamFile>&& bamFiles);
    PbiFilterCompositeBamReader(const PbiFilter& filter, const DataSet& dataset);

public:
    bool GetNext(BamRecord& record);

public:
    PbiFilterCompositeBamReader& Filter(const PbiFilter& filter);

private:
    void UpdateSort(void);

private:
    container_type mergeQueue_;
    std::vector<std::string> filenames_;
};

class SequentialCompositeBamReader
{
public:
    SequentialCompositeBamReader(const std::vector<BamFile>& bamFiles);
    SequentialCompositeBamReader(std::vector<BamFile>&& bamFiles);
    SequentialCompositeBamReader(const DataSet& dataset);

public:
    bool GetNext(BamRecord& record);

private:
    std::deque<std::unique_ptr<BamReader> > readers_;
};

} // namespace internal
} // namespace BAM
} // namespace PacBio

#include "CompositeBamReader.inl"

#endif // COMPOSITEBAMREADER_H
