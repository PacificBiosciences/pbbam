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
/// \file PbiBasicTypes.h
/// \brief Defines the basic data structures used in PBI lookups.
//
// Author: Derek Barnett

#ifndef PBIBASICTYPES_H
#define PBIBASICTYPES_H

#include "pbbam/Compare.h"
#include "pbbam/Config.h"
#include <cstddef>
#include <cstdint>
#include <deque>
#include <utility>
#include <vector>

namespace PacBio {
namespace BAM {

/// \brief The IndexResultBlock class represents a contiguous group of records
///        returned from a PBI lookup.
///
/// Contiguous reads that satisfy a PBI lookup query will be merged down into a
/// single block. This helps to minimize the number of seeks in subsequent read
/// operations.
///
/// An PBI-enabled reader or query can iterate over a list of IndexResultBlocks;
/// for each block, seeking to the first record and then sequentially reading
/// 'numReads' consecutive records before needing to seek again.
///
struct PBBAM_EXPORT IndexResultBlock
{
public:
    IndexResultBlock() = default;
    IndexResultBlock(size_t idx, size_t numReads);

public:
    bool operator==(const IndexResultBlock& other) const;
    bool operator!=(const IndexResultBlock& other) const;

public:
    size_t  firstIndex_ = 0;     ///< index of block's first record in BAM/PBI files (e.g. i-th record)
    size_t  numReads_ = 0;       ///< number of reads in this block
    int64_t virtualOffset_ = -1;  ///< virtual offset of first record in this block
};

/// \brief container of PBI result blocks
///
using IndexResultBlocks = std::deque<IndexResultBlock>;

/// \brief container of raw PBI indices
///
/// This is the primary result of PbiFilter -associated classes. This raw list
/// can participate in set operations (union, intersect) for compound filters,
/// and then be merged down into IndexResultBlocks for actual data file
/// random-access.
///
using IndexList = std::vector<size_t>;

/// \brief pair representing a range of PBI indices: where interval
///        is [first, second)
///
/// Used primarily by the PBI's CoordinateSortedData components.
///
/// \sa PbiReferenceEntry, PbiRawReferenceData, & ReferenceLookupData
///
using IndexRange = std::pair<size_t, size_t>;

} // namespace BAM
} // namespace PacBio

#include "pbbam/internal/PbiBasicTypes.inl"

#endif // PBIBASICTYPES_H
