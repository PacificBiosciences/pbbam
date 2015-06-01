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

#ifndef BAMRECORDSORT_H
#define BAMRECORDSORT_H

#include "pbbam/BamRecord.h"
#include <functional>
#include <cassert>

namespace PacBio {
namespace BAM {
namespace internal {

enum class SortOrder {
    Ascending = 0
  , Descending
};

template<typename ElemType>
inline bool sort_helper(const SortOrder& order,
                        const ElemType& lhs,
                        const ElemType& rhs)
{
    switch ( order ) {
        case SortOrder::Ascending   : { std::less<ElemType> comp;    return comp(lhs, rhs); }
        case SortOrder::Descending  : { std::greater<ElemType> comp; return comp(lhs, rhs); }
        default :
            assert(false);
    }
    return false; // <-- unreachable
}

typedef std::binary_function<BamRecord, BamRecord, bool> BamRecordSortBase;

struct Unsorted : public BamRecordSortBase
{
public:
    Unsorted(const SortOrder& order = SortOrder::Ascending)
    { (void)order; }

    bool operator()(const BamRecord& lhs, const BamRecord& rhs)
    { (void)lhs; (void)rhs; return false; }
};

struct ByQName : public BamRecordSortBase
{
public:
    ByQName(const SortOrder& order = SortOrder::Ascending)
        : m_order(order)
    { }

    bool operator()(const BamRecord& lhs, const BamRecord& rhs)
    { return sort_helper(m_order, lhs.FullName(), rhs.FullName()); }

private:
    const SortOrder m_order;
};

struct ByPosition : public BamRecordSortBase
{
public:
    ByPosition(const SortOrder& order = SortOrder::Ascending)
        : m_order(order)
    { }

    // comparison function
    bool operator()(const BamRecord& lhs, const BamRecord& rhs) {

        const int32_t lhsId = lhs.ReferenceId();
        const int32_t rhsId = rhs.ReferenceId();

        // force unmapped aligmnents to end
        if ( lhsId == -1 ) return false;
        if ( rhsId == -1 ) return true;

        // if on same reference, sort on position
        if ( lhsId == rhsId )
            return sort_helper(m_order, lhs.ReferenceStart(), rhs.ReferenceStart());

        // otherwise sort on reference ID
        return sort_helper(m_order, lhsId, rhsId);
    }

private:
    const SortOrder m_order;
};

struct ByZmw : public BamRecordSortBase {
public:
    ByZmw(const SortOrder& order = SortOrder::Ascending) : m_order(order) { }

    bool operator()(const BamRecord& lhs, const BamRecord& rhs)
    { return sort_helper(m_order, lhs.HoleNumber(), rhs.HoleNumber()); }

private:
    const SortOrder m_order;
};

} // namespace internal
} // namespace BAM
} // namespace PacBio

#endif // BAMRECORDSORT_H
