// File Description
/// \file PbiBasicTypes.h
/// \brief Defines the basic data structures used in PBI lookups.
//
// Author: Derek Barnett

#ifndef PBIBASICTYPES_H
#define PBIBASICTYPES_H

#include <cstddef>
#include <cstdint>
#include <deque>
#include <utility>
#include <vector>
#include "pbbam/Compare.h"
#include "pbbam/Config.h"

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
    IndexResultBlock(size_t idx, size_t numReads);

    IndexResultBlock() = default;
    IndexResultBlock(const IndexResultBlock&) = default;
    IndexResultBlock(IndexResultBlock&&) = default;
    IndexResultBlock& operator=(const IndexResultBlock&) = default;
    IndexResultBlock& operator=(IndexResultBlock&&) = default;
    ~IndexResultBlock() = default;

public:
    bool operator==(const IndexResultBlock& other) const;
    bool operator!=(const IndexResultBlock& other) const;

public:
    size_t firstIndex_ = 0;  ///< index of block's first record in BAM/PBI files (e.g. i-th record)
    size_t numReads_ = 0;    ///< number of reads in this block
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

}  // namespace BAM
}  // namespace PacBio

#include "pbbam/internal/PbiBasicTypes.inl"

#endif  // PBIBASICTYPES_H
