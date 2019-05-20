// File Description
/// \file CompositeBamReader.inl
/// \brief Inline implementations for the composite BAM readers, for
///        working with multiple input files.
//
// Author: Derek Barnett

#include "pbbam/CompositeBamReader.h"

#include <algorithm>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>

namespace PacBio {
namespace BAM {
namespace internal {

// -----------------------------------
// Merging helpers
// -----------------------------------

inline CompositeMergeItem::CompositeMergeItem(std::unique_ptr<BamReader> rdr)
    : reader{std::move(rdr)}
{
}

inline CompositeMergeItem::CompositeMergeItem(std::unique_ptr<BamReader> rdr, BamRecord rec)
    : reader{std::move(rdr)}, record{std::move(rec)}
{
}

template <typename CompareType>
inline bool CompositeMergeItemSorter<CompareType>::operator()(const CompositeMergeItem& lhs,
                                                              const CompositeMergeItem& rhs) const
{
    const auto& l = lhs.record;
    const auto& r = rhs.record;
    return CompareType()(l, r);
}

}  // namespace internal

// -----------------------------------
// GenomicIntervalCompositeBamReader
// -----------------------------------

inline GenomicIntervalCompositeBamReader::GenomicIntervalCompositeBamReader(
    const std::vector<BamFile>& bamFiles)
    : GenomicIntervalCompositeBamReader{bamFiles, MakeBaiIndexCache(bamFiles)}
{
}

inline GenomicIntervalCompositeBamReader::GenomicIntervalCompositeBamReader(
    const std::vector<BamFile>& bamFiles, const BaiIndexCache& cache)
{
    indexCache_ = cache;

    filenames_.reserve(bamFiles.size());
    for (const auto& bamFile : bamFiles)
        filenames_.push_back(bamFile.Filename());
}

inline GenomicIntervalCompositeBamReader::GenomicIntervalCompositeBamReader(const DataSet& dataset)
    : GenomicIntervalCompositeBamReader{dataset.BamFiles()}
{
}

inline GenomicIntervalCompositeBamReader::GenomicIntervalCompositeBamReader(
    const DataSet& dataset, const BaiIndexCache& cache)
    : GenomicIntervalCompositeBamReader{dataset.BamFiles(), cache}
{
}

inline GenomicIntervalCompositeBamReader::GenomicIntervalCompositeBamReader(
    const GenomicInterval& interval, const std::vector<BamFile>& bamFiles)
    : GenomicIntervalCompositeBamReader{interval, bamFiles, MakeBaiIndexCache(bamFiles)}
{
}

inline GenomicIntervalCompositeBamReader::GenomicIntervalCompositeBamReader(
    const GenomicInterval& interval, const std::vector<BamFile>& bamFiles,
    const BaiIndexCache& cache)
    : GenomicIntervalCompositeBamReader{bamFiles, cache}
{
    Interval(interval);
}

inline GenomicIntervalCompositeBamReader::GenomicIntervalCompositeBamReader(
    const GenomicInterval& interval, const DataSet& dataset)
    : GenomicIntervalCompositeBamReader{interval, dataset.BamFiles()}
{
}

inline GenomicIntervalCompositeBamReader::GenomicIntervalCompositeBamReader(
    const GenomicInterval& interval, const DataSet& dataset, const BaiIndexCache& cache)
    : GenomicIntervalCompositeBamReader{interval, dataset.BamFiles(), cache}
{
}

inline bool GenomicIntervalCompositeBamReader::GetNext(BamRecord& record)
{
    // nothing left to read
    if (mergeItems_.empty()) return false;

    // non-destructive 'pop' of first item from queue
    auto firstIter = mergeItems_.begin();
    auto firstItem =
        internal::CompositeMergeItem{std::move(firstIter->reader), std::move(firstIter->record)};
    mergeItems_.pop_front();

    // store its record in our output record
    std::swap(record, firstItem.record);

    // try fetch 'next' from first item's reader
    // if successful, re-insert it into container & re-sort on our new values
    // otherwise, this item will go out of scope & reader destroyed
    if (firstItem.reader->GetNext(firstItem.record)) {
        mergeItems_.push_front(std::move(firstItem));
        UpdateSort();
    }

    // return success
    return true;
}

inline const GenomicInterval& GenomicIntervalCompositeBamReader::Interval() const
{
    return interval_;
}

inline GenomicIntervalCompositeBamReader& GenomicIntervalCompositeBamReader::Interval(
    const GenomicInterval& interval)
{
    // reset readers
    mergeItems_.clear();

    // create readers for files
    std::deque<internal::CompositeMergeItem> updatedMergeItems;
    std::vector<std::string> missingBai;
    for (size_t i = 0; i < filenames_.size(); ++i) {
        const BamFile bamFile{filenames_.at(i)};
        if (bamFile.StandardIndexExists()) {
            internal::CompositeMergeItem item{std::unique_ptr<BamReader>{
                new BaiIndexedBamReader{interval, std::move(bamFile), indexCache_->at(i)}}};
            if (item.reader->GetNext(item.record)) updatedMergeItems.push_back(std::move(item));
            // else not an error, simply no data matching interval
        } else {
            // maybe handle PBI-backed interval searches if BAI missing, but for now treat as error
            missingBai.push_back(bamFile.Filename());
        }
    }

    // throw if any files missing BAI
    if (!missingBai.empty()) {
        std::ostringstream e;
        e << "GenomicIntervalCompositeBamReader: failed to open because the following files are "
             "missing a *.bai index:\n";
        for (const auto& fn : missingBai)
            e << "  " << fn << '\n';
        throw std::runtime_error{e.str()};
    }

    // update our actual container and return
    mergeItems_ = std::move(updatedMergeItems);
    UpdateSort();
    return *this;
}

struct OrderByPosition
{
    static inline bool less_than(const BamRecord& lhs, const BamRecord& rhs)
    {
        const int32_t lhsId = lhs.ReferenceId();
        const int32_t rhsId = rhs.ReferenceId();
        if (lhsId == -1) return false;
        if (rhsId == -1) return true;

        if (lhsId == rhsId)
            return lhs.ReferenceStart() < rhs.ReferenceStart();
        else
            return lhsId < rhsId;
    }

    static inline bool equals(const BamRecord& lhs, const BamRecord& rhs)
    {
        return lhs.ReferenceId() == rhs.ReferenceId() &&
               lhs.ReferenceStart() == rhs.ReferenceStart();
    }
};

struct PositionSorter
    : std::binary_function<internal::CompositeMergeItem, internal::CompositeMergeItem, bool>
{
    bool operator()(const internal::CompositeMergeItem& lhs,
                    const internal::CompositeMergeItem& rhs) const
    {
        const BamRecord& l = lhs.record;
        const BamRecord& r = rhs.record;
        return OrderByPosition::less_than(l, r);
    }
};

inline void GenomicIntervalCompositeBamReader::UpdateSort()
{
    std::sort(mergeItems_.begin(), mergeItems_.end(), PositionSorter{});
}

// ------------------------------
// PbiRequestCompositeBamReader
// ------------------------------

template <typename OrderByType>
inline PbiFilterCompositeBamReader<OrderByType>::PbiFilterCompositeBamReader(
    const PbiFilter& filter, const std::vector<BamFile>& bamFiles)
    : PbiFilterCompositeBamReader{filter, bamFiles, MakePbiIndexCache(bamFiles)}
{
}

template <typename OrderByType>
inline PbiFilterCompositeBamReader<OrderByType>::PbiFilterCompositeBamReader(
    const PbiFilter& filter, const std::vector<BamFile>& bamFiles, const PbiIndexCache& cache)
    : indexCache_{cache}, numReads_{0}
{
    filenames_.reserve(bamFiles.size());
    for (const auto& bamFile : bamFiles)
        filenames_.push_back(bamFile.Filename());
    Filter(filter);
}

template <typename OrderByType>
inline PbiFilterCompositeBamReader<OrderByType>::PbiFilterCompositeBamReader(
    const PbiFilter& filter, const DataSet& dataset)
    : PbiFilterCompositeBamReader{filter, dataset.BamFiles()}
{
}

template <typename OrderByType>
inline PbiFilterCompositeBamReader<OrderByType>::PbiFilterCompositeBamReader(
    const PbiFilter& filter, const DataSet& dataset, const PbiIndexCache& cache)
    : PbiFilterCompositeBamReader{filter, dataset.BamFiles(), cache}
{
}

template <typename OrderByType>
inline bool PbiFilterCompositeBamReader<OrderByType>::GetNext(BamRecord& record)
{
    // nothing left to read
    if (mergeQueue_.empty()) return false;

    // non-destructive 'pop' of first item from queue
    auto firstIter = mergeQueue_.begin();
    value_type firstItem{std::move(firstIter->reader), std::move(firstIter->record)};
    mergeQueue_.pop_front();

    // store its record in our output record
    std::swap(record, firstItem.record);

    // try fetch 'next' from first item's reader
    // if successful, re-insert it into container & re-sort on our new values
    // otherwise, this item will go out of scope & reader destroyed
    if (firstItem.reader->GetNext(firstItem.record)) {
        mergeQueue_.push_front(std::move(firstItem));
        UpdateSort();
    }

    // return success
    return true;
}

template <typename OrderByType>
inline PbiFilterCompositeBamReader<OrderByType>& PbiFilterCompositeBamReader<OrderByType>::Filter(
    const PbiFilter& filter)
{
    // std::cerr << "PbiFilterCompositeBamReader<OrderByType>::Filter()\n";

    // reset reader queue
    mergeQueue_.clear();

    // create readers for files
    container_type updatedMergeItems;
    std::vector<std::string> missingPbi;
    for (size_t i = 0; i < filenames_.size(); ++i) {
        const BamFile bamFile{filenames_.at(i)};
        if (bamFile.PacBioIndexExists()) {
            auto item = internal::CompositeMergeItem{std::unique_ptr<BamReader>{
                new PbiIndexedBamReader{filter, std::move(bamFile), indexCache_->at(i)}}};
            if (item.reader->GetNext(item.record)) updatedMergeItems.push_back(std::move(item));
            // else not an error, simply no data matching filter
        } else
            missingPbi.push_back(filenames_.at(i));
    }

    // throw if any files missing PBI
    if (!missingPbi.empty()) {
        std::ostringstream e;
        e << "PbiFilterCompositeBamReader: failed to open because the following files are "
             "missing a *.pbi index:\n";
        for (const auto& fn : missingPbi)
            e << "  " << fn << '\n';
        throw std::runtime_error{e.str()};
    }

    // update our actual container, store num matching reads, sort & and return
    mergeQueue_ = std::move(updatedMergeItems);
    numReads_ = 0;
    for (const auto& item : mergeQueue_) {
        auto* pbiReader = dynamic_cast<PbiIndexedBamReader*>(item.reader.get());
        numReads_ += pbiReader->NumReads();
    }
    UpdateSort();
    return *this;
}

template <typename OrderByType>
inline uint32_t PbiFilterCompositeBamReader<OrderByType>::NumReads() const
{
    return numReads_;
}

template <typename OrderByType>
inline void PbiFilterCompositeBamReader<OrderByType>::UpdateSort()
{
    std::stable_sort(mergeQueue_.begin(), mergeQueue_.end(), merge_sorter_type{});
}

// ------------------------------
// SequentialCompositeBamReader
// ------------------------------

inline SequentialCompositeBamReader::SequentialCompositeBamReader(std::vector<BamFile> bamFiles)
{
    for (auto&& bamFile : bamFiles)
        readers_.emplace_back(std::make_unique<BamReader>(std::move(bamFile)));
}

inline SequentialCompositeBamReader::SequentialCompositeBamReader(const DataSet& dataset)
    : SequentialCompositeBamReader{dataset.BamFiles()}
{
}

inline bool SequentialCompositeBamReader::GetNext(BamRecord& record)
{
    // try first reader, if successful return true
    // else pop reader and try next, until all readers exhausted
    while (!readers_.empty()) {
        auto& reader = readers_.front();
        if (reader->GetNext(record))
            return true;
        else
            readers_.pop_front();
    }

    // no readers available
    return false;
}

}  // namespace BAM
}  // namespace PacBio
