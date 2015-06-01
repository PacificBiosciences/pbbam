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

#ifndef MERGESTRATEGY_H
#define MERGESTRATEGY_H

#include "pbbam/BamRecord.h"
#include "pbbam/internal/IMergeStrategy.h"
#include "pbbam/internal/MergeItem.h"
#include <functional>
#include <set>
#include <vector>
#include <cassert>

namespace PacBio {
namespace BAM {
namespace internal {

template<typename Compare>
struct MergeItemSorter : public std::binary_function<MergeItem, MergeItem, bool>
{
public:
    MergeItemSorter(const Compare& comp = Compare())
        : comp_(comp)
    { }

    bool operator()(const MergeItem& lhs, const MergeItem& rhs) {
        const BamRecord& l = lhs.record_;
        const BamRecord& r = rhs.record_;
        return comp_(l, r);
    }

private:
    Compare comp_;
};

template<typename Compare>
struct GroupMergeItemSorter : public std::binary_function<GroupMergeItem, GroupMergeItem, bool>
{
public:
    GroupMergeItemSorter(const Compare& comp = Compare())
        : comp_(comp)
    { }

    bool operator()(const GroupMergeItem& lhs, const GroupMergeItem& rhs) {
        if ( lhs.record_.empty())
            return false;
        if ( rhs.record_.empty())
            return true;
        assert(!lhs.record_.empty());
        assert(!rhs.record_.empty());
        const BamRecord& l = lhs.record_.front();
        const BamRecord& r = rhs.record_.front();
        return comp_(l, r);
    }

private:
    Compare comp_;
};

template<typename Compare>
class MergeStrategy : public IMergeStrategy
{
public:
    MergeStrategy(const std::vector<FileIterPtr>& iters);
    bool GetNext(BamRecord& record);
private:
    std::multiset<MergeItem, MergeItemSorter<Compare> > mergeItems_;
};

template<typename Compare>
class GroupMergeStrategy : public IGroupMergeStrategy
{
public:
    GroupMergeStrategy(const std::vector<FileIterPtr>& iters);
    bool GetNext(std::vector<BamRecord>& records);
private:
    GroupMergeItem nextItem_;
    std::multiset<GroupMergeItem, GroupMergeItemSorter<Compare> > mergeItems_;
};

// -----------------------
// MergeStrategy
// -----------------------

template<typename Compare>
inline MergeStrategy<Compare>::MergeStrategy(const std::vector<FileIterPtr>& iters)
    : IMergeStrategy()
{
    BamRecord r;
    for (FileIterPtr iter : iters) {
        if (iter->GetNext(r)) {
            MergeItem item(r, iter);
            mergeItems_.insert(item);
        }
    }
}

template<typename Compare>
inline bool MergeStrategy<Compare>::GetNext(BamRecord& record)
{
    if (mergeItems_.empty())
        return false;

    // pop first merge item & record
    auto firstIter = mergeItems_.begin();
    MergeItem firstItem = (*firstIter);
    mergeItems_.erase(firstIter);
    record = firstItem.record_;

    // try fetch iter's next (if failed, do not replace)
    if (firstItem.iter_->GetNext(firstItem.record_))
        mergeItems_.insert(firstItem);
    return true;
}

// -----------------------
// GroupMergeStrategy
// -----------------------

template<typename Compare>
inline GroupMergeStrategy<Compare>::GroupMergeStrategy(const std::vector<FileIterPtr>& iters)
    : IGroupMergeStrategy()
{
    std::vector<BamRecord> r;
    for (FileIterPtr iter : iters) {
        if (iter->GetNext(r)) {
            GroupMergeItem item(r, iter);
            mergeItems_.insert(item);
        }
    }
    if (!mergeItems_.empty()) {
        auto firstIter = mergeItems_.begin();
        nextItem_ = (*firstIter);
        mergeItems_.erase(firstIter);
    }
}

template<typename Compare>
inline bool GroupMergeStrategy<Compare>::GetNext(std::vector<BamRecord>& records)
{
    records.clear();
    if (nextItem_.IsNull())
        return false;

    // append "nextItem" records
    records = nextItem_.record_;

    // try fetch iter's next (if failed, do not replace)
    if (nextItem_.iter_->GetNext(nextItem_.record_))
        mergeItems_.insert(nextItem_);
    else
        nextItem_ = GroupMergeItem();

    while (!mergeItems_.empty()) {

        // pop first merge item
        auto firstIter = mergeItems_.begin();
        GroupMergeItem firstItem = (*firstIter);
        mergeItems_.erase(firstIter);

        // if first item has records
        if (!firstItem.record_.empty()) {

            // if first block to store
            if (records.empty())
                records = firstItem.record_;

            // else see if we match current group
            else {
                const BamRecord& lhs = records.front();
                const BamRecord& rhs = firstItem.record_.front();

                // if match, append to output & fetch next
                if (firstItem.iter_->InSameGroup(lhs, rhs)) {

                    for (const BamRecord& r : firstItem.record_)
                        records.push_back(r);
                    if (firstItem.iter_->GetNext(firstItem.record_))
                        mergeItems_.insert(firstItem);
                }

                // no match, item becomes the "next item" to use
                else {
                    nextItem_ = firstItem;
                    break;
                }
            }
        }

        // first item has no records, try fetch next
        else {
            if (firstItem.iter_->GetNext(firstItem.record_))
                mergeItems_.insert(firstItem);
        }
    }

    return true;
}

} // namespace internal
} // namespace BAM
} // namespace PacBio

#endif // MERGESTRATEGY_H
