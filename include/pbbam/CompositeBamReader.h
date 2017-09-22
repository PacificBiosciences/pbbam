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
/// \file CompositeBamReader.h
/// \brief Defines the composite BAM readers, for working with multiple input
///       files.
//
// Author: Derek Barnett

#ifndef COMPOSITEBAMREADER_H
#define COMPOSITEBAMREADER_H

#include "pbbam/BaiIndexedBamReader.h"
#include "pbbam/BamFile.h"
#include "pbbam/BamHeader.h"
#include "pbbam/BamReader.h"
#include "pbbam/BamRecord.h"
#include "pbbam/Config.h"
#include "pbbam/DataSet.h"
#include "pbbam/GenomicInterval.h"
#include "pbbam/PbiIndexedBamReader.h"
#include <deque>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

namespace internal {

/// \internal
/// \brief The CompositeMergeItem class provides a helper struct for composite
///        readers, containing a single-file reader and its "next" record.
///
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
    ~CompositeMergeItem() = default;
};

/// \internal
/// \brief The CompositeMergeItemSorter class provides a helper function object
///        for ordering composite reader results.
///
/// Essentially just exracts a BamRecord from its parent CompositeMergeItem for
/// further checks.
///
template<typename CompareType>
struct CompositeMergeItemSorter : public std::function<bool(const CompositeMergeItem&,
                                                            const CompositeMergeItem&)>
{
    bool operator()(const CompositeMergeItem& lhs,
                    const CompositeMergeItem& rhs);
};

} // namespace internal

/// \brief The GenomicIntervalCompositeBamReader class provides read access to
///        multipe %BAM files, limiting results to a genomic region.
///
/// Requires a ".bai" file for each input %BAM file.
///
/// Results will be returned in order of genomic coordinate (first by reference
/// ID, then by position).
///
class PBBAM_EXPORT GenomicIntervalCompositeBamReader
{
public:
    /// \name Contstructors & Related Methods
    /// \{

    GenomicIntervalCompositeBamReader(const GenomicInterval& interval,
                                      const std::vector<BamFile>& bamFiles);
    GenomicIntervalCompositeBamReader(const GenomicInterval& interval,
                                      std::vector<BamFile>&& bamFiles);
    GenomicIntervalCompositeBamReader(const GenomicInterval& interval,
                                      const DataSet& dataset);

    /// \}

public:
    /// \name Data Access
    /// \{

    /// Fetches next BAM record in the interval specified, storing in \p record
    ///
    /// \param[out] record
    /// \returns true on success, false if no more data available.
    ///
    bool GetNext(BamRecord& record);

    /// Sets a new genomic interval of interest.
    ///
    /// \returns reference to this reader
    ///
    GenomicIntervalCompositeBamReader& Interval(const GenomicInterval& interval);

    /// \returns the current specified interval
    ///
    const GenomicInterval& Interval() const;

    /// \}

private:
    void UpdateSort();

private:
    GenomicInterval interval_;
    std::deque<internal::CompositeMergeItem> mergeItems_;
    std::vector<std::string> filenames_;
};

/// \brief Provides read access to multipe %BAM files, limiting results to those
///        passing a PbiFilter.
///
/// Requires a ".pbi" file for each input %BAM file.
///
/// \note The template parameter OrderByType is not fully implemented at this
///       time. Use of comparison functor (e.g. Compare::Zmw) for this will
///       currently result in the proper "next" value <b> at each iteration
///       step, independently, but not over the full data set. </b> If all
///       files' "order-by" data values are accessible in increasing order
///       within each file, then the expected ordering will be observed,
///       However, if these data are not sorted within a file, the final results
///       will appear unordered. \n
///       \n
///           Example:\n
///           file 1: { 1, 5, 2, 6 } \n
///           file 2: { 3, 8, 4, 7 } \n
///           results: { 1, 3, 5, 2, 6, 8, 4, 7 } \n
///       \n
///       This a known issue and will be addressed in a future update. But in
///       the meantime, use of Compare::None as the OrderByType is recommended,
///       to explicitly indicate that no particular ordering is expected.
///
template<typename OrderByType>
class PBBAM_EXPORT PbiFilterCompositeBamReader
{
public:
    using value_type        = internal::CompositeMergeItem;
    using merge_sorter_type = internal::CompositeMergeItemSorter<OrderByType>;
    using container_type    = std::deque<value_type>;
    using iterator          = typename container_type::iterator;
    using const_iterator    = typename container_type::const_iterator;

public:
    /// \name Contstructors & Related Methods
    /// \{

    PbiFilterCompositeBamReader(const PbiFilter& filter,
                                const std::vector<BamFile>& bamFiles);
    PbiFilterCompositeBamReader(const PbiFilter& filter,
                                std::vector<BamFile>&& bamFiles);
    PbiFilterCompositeBamReader(const PbiFilter& filter,
                                const DataSet& dataset);

    /// \}

public:
    /// \name Data Access
    /// \{

    /// Fetches next BAM record in the interval specified.
    ///
    /// \returns true on success, false if no more data available.
    ///
    bool GetNext(BamRecord& record);

    /// Sets a new PBI filter
    ///
    /// \returns reference to this reader
    ///
    PbiFilterCompositeBamReader& Filter(const PbiFilter& filter);

    uint32_t NumReads() const;

    /// \}

private:
    void UpdateSort();

private:
    container_type mergeQueue_;
    std::vector<std::string> filenames_;
    uint32_t numReads_;
};

/// \brief The SequentialCompositeBamReader class provides read access to
///        multiple %BAM files, reading through the entire contents of each
///        file.
///
/// Input files will be accessed in the order provided to the constructor. Each
/// file's contents will be exhausted before moving on to the next one (as
/// opposed to a "round-robin" scheme).
///
class PBBAM_EXPORT SequentialCompositeBamReader
{
public:
    /// \name Contstructors & Related Methods
    /// \{

    SequentialCompositeBamReader(const std::vector<BamFile>& bamFiles);
    SequentialCompositeBamReader(std::vector<BamFile>&& bamFiles);
    SequentialCompositeBamReader(const DataSet& dataset);

    /// \}

public:
    /// \name Data Access
    /// \{

    /// Fetches next BAM record from the .
    ///
    /// \returns true on success, false if no more data available.
    ///
    bool GetNext(BamRecord& record);

    /// \}

private:
    std::deque<std::unique_ptr<BamReader> > readers_;
};

} // namespace BAM
} // namespace PacBio

#include "pbbam/internal/CompositeBamReader.inl"

#endif // COMPOSITEBAMREADER_H
