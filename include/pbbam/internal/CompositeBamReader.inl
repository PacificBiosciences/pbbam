// File Description
/// \file CompositeBamReader.inl
/// \brief Inline implementations for the composite BAM readers, for
///        working with multiple input files.
//
// Author: Derek Barnett

#include "pbbam/CompositeBamReader.h"
#include <algorithm>
#include <set>
#include <sstream>
#include <stdexcept>

#include "pbbam/MakeUnique.h"

namespace PacBio {
namespace BAM {
namespace internal {

// -----------------------------------
// Merging helpers
// -----------------------------------

inline CompositeMergeItem::CompositeMergeItem(std::unique_ptr<BamReader> rdr)
    : reader{std::move(rdr)}
{ }

inline CompositeMergeItem::CompositeMergeItem(std::unique_ptr<BamReader> rdr,
                                              BamRecord rec)
    : reader{std::move(rdr)}
    , record{std::move(rec)}
{ }

template<typename CompareType>
inline bool CompositeMergeItemSorter<CompareType>::operator()(const CompositeMergeItem& lhs,
                                                              const CompositeMergeItem& rhs)
{
    const auto& l = lhs.record;
    const auto& r = rhs.record;
    return CompareType()(l, r);
}

} // namespace internal

// -----------------------------------
// GenomicIntervalCompositeBamReader
// -----------------------------------

inline GenomicIntervalCompositeBamReader::GenomicIntervalCompositeBamReader(const GenomicInterval& interval,
                                                                            const std::vector<BamFile>& bamFiles)
{
    filenames_.reserve(bamFiles.size());
    for(const auto& bamFile : bamFiles)
        filenames_.push_back(bamFile.Filename());
    Interval(interval);
}

inline GenomicIntervalCompositeBamReader::GenomicIntervalCompositeBamReader(const GenomicInterval& interval,
                                                                            const DataSet& dataset)
    : GenomicIntervalCompositeBamReader{interval, dataset.BamFiles()}
{ }

inline bool GenomicIntervalCompositeBamReader::GetNext(BamRecord& record)
{
    // nothing left to read
    if (mergeItems_.empty())
        return false;

    // non-destructive 'pop' of first item from queue
    auto firstIter = mergeItems_.begin();
    auto firstItem = internal::CompositeMergeItem{ std::move(firstIter->reader), std::move(firstIter->record) };
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
{ return interval_; }

inline GenomicIntervalCompositeBamReader& GenomicIntervalCompositeBamReader::Interval(const GenomicInterval& interval)
{
    std::deque<internal::CompositeMergeItem> updatedMergeItems;
    std::set<std::string> filesToCreate{filenames_.cbegin(), filenames_.cend()};

    // update existing readers
    while (!mergeItems_.empty()) {

        // non-destructive 'pop' of first item from queue
        auto firstIter = mergeItems_.begin();
        internal::CompositeMergeItem firstItem{ std::move(firstIter->reader), std::move(firstIter->record) };
        mergeItems_.pop_front();

        // reset interval
        auto* baiReader = dynamic_cast<BaiIndexedBamReader*>(firstItem.reader.get());
        assert(baiReader);
        baiReader->Interval(interval);

        // try fetch 'next' from first item's reader
        // if successful, re-insert it into container & re-sort on our new values
        // otherwise, this item will go out of scope & reader destroyed
        if (firstItem.reader->GetNext(firstItem.record)) {
            updatedMergeItems.push_front(std::move(firstItem));
            filesToCreate.erase(firstItem.reader->Filename());
        }
    }

    // create readers for files that were not 'active' for the previous
    std::vector<std::string> missingBai;
    for (auto&& fn : filesToCreate) {
        BamFile bamFile{ fn };
        if (bamFile.StandardIndexExists()) {
            internal::CompositeMergeItem item{ std::unique_ptr<BamReader>{ new BaiIndexedBamReader{ interval, std::move(bamFile) } } };
            if (item.reader->GetNext(item.record))
                updatedMergeItems.push_back(std::move(item));
            // else not an error, simply no data matching interval
        }
        else {
            // maybe handle PBI-backed interval searches if BAI missing, but for now treat as error
            missingBai.push_back(bamFile.Filename());
        }
    }

    // throw if any files missing BAI
    if (!missingBai.empty()) {
        std::ostringstream e;
        e << "failed to open GenomicIntervalCompositeBamReader because the following files are missing a BAI file:\n";
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
        else return lhsId < rhsId;
    }

    static inline bool equals(const BamRecord& lhs, const BamRecord& rhs)
    {
        return lhs.ReferenceId() == rhs.ReferenceId() &&
               lhs.ReferenceStart() == rhs.ReferenceStart();
    }
};

struct PositionSorter : std::binary_function<internal::CompositeMergeItem, internal::CompositeMergeItem, bool>
{
    bool operator()(const internal::CompositeMergeItem& lhs,
                    const internal::CompositeMergeItem& rhs)
    {
        const BamRecord& l = lhs.record;
        const BamRecord& r = rhs.record;
        return OrderByPosition::less_than(l, r);
    }
};

inline void GenomicIntervalCompositeBamReader::UpdateSort()
{ std::sort(mergeItems_.begin(), mergeItems_.end(), PositionSorter{ }); }

// ------------------------------
// PbiRequestCompositeBamReader
// ------------------------------

template<typename OrderByType>
inline PbiFilterCompositeBamReader<OrderByType>::PbiFilterCompositeBamReader(const PbiFilter& filter,
                                                                             const std::vector<BamFile>& bamFiles)
    : numReads_{0}
{
    filenames_.reserve(bamFiles.size());
    for(const auto& bamFile : bamFiles)
        filenames_.push_back(bamFile.Filename());
    Filter(filter);
}

template<typename OrderByType>
inline PbiFilterCompositeBamReader<OrderByType>::PbiFilterCompositeBamReader(const PbiFilter& filter,
                                                                             const DataSet& dataset)
    : PbiFilterCompositeBamReader{filter, dataset.BamFiles()}
{ }

template<typename OrderByType>
inline bool PbiFilterCompositeBamReader<OrderByType>::GetNext(BamRecord& record)
{
    // nothing left to read
    if (mergeQueue_.empty())
        return false;

    // non-destructive 'pop' of first item from queue
    auto firstIter = mergeQueue_.begin();
    value_type firstItem{ std::move(firstIter->reader), std::move(firstIter->record) };
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

template<typename OrderByType>
inline PbiFilterCompositeBamReader<OrderByType>&
PbiFilterCompositeBamReader<OrderByType>::Filter(const PbiFilter& filter)
{
    container_type updatedMergeItems;
    std::set<std::string> filesToCreate{ filenames_.cbegin(), filenames_.cend() };

    // update existing readers
    while (!mergeQueue_.empty()) {

        // non-destructive 'pop' of first item from queue
        auto firstIter = mergeQueue_.begin();
        internal::CompositeMergeItem firstItem{ std::move(firstIter->reader), std::move(firstIter->record) };
        mergeQueue_.pop_front();

        // reset request
        auto* pbiReader = dynamic_cast<PbiIndexedBamReader*>(firstItem.reader.get());
        assert(pbiReader);
        pbiReader->Filter(filter);

        // try fetch 'next' from first item's reader
        // if successful, re-insert it into container & re-sort on our new values
        // otherwise, this item will go out of scope & reader destroyed
        if (firstItem.reader->GetNext(firstItem.record)) {
            updatedMergeItems.push_front(std::move(firstItem));
            filesToCreate.erase(firstItem.reader->Filename());
        }
    }

    // create readers for files that were not 'active' for the previous
    std::vector<std::string> missingPbi;
    for (auto&& fn : filesToCreate) {
        const BamFile bamFile{ fn };
        if (bamFile.PacBioIndexExists()) {
            auto item = internal::CompositeMergeItem{ std::unique_ptr<BamReader>{ new PbiIndexedBamReader{ filter, std::move(bamFile) } } };
            if (item.reader->GetNext(item.record))
                updatedMergeItems.push_back(std::move(item));
            // else not an error, simply no data matching filter
        }
        else
            missingPbi.push_back(fn);
    }

    // throw if any files missing PBI
    if (!missingPbi.empty()) {
        std::ostringstream e;
        e << "failed to open PbiFilterCompositeBamReader because the following files are missing a PBI file:\n";
        for (const auto& fn : missingPbi)
            e << "  " << fn << '\n';
        throw std::runtime_error{e.str()};
    }


    // update our actual container, store num matching reads, sort & and return
    mergeQueue_ = std::move(updatedMergeItems);

    numReads_ = 0;
    for (const auto& item : mergeQueue_)
    {
        auto* pbiReader = dynamic_cast<PbiIndexedBamReader*>(item.reader.get());
        numReads_ += pbiReader->NumReads();
    }

    UpdateSort();
    return *this;
}

template<typename OrderByType>
inline uint32_t PbiFilterCompositeBamReader<OrderByType>::NumReads() const
{
    return numReads_;
}

template<typename OrderByType>
inline void PbiFilterCompositeBamReader<OrderByType>::UpdateSort()
{ std::stable_sort(mergeQueue_.begin(), mergeQueue_.end(), merge_sorter_type{}); }

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
{ }

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

} // namespace BAM
} // namespace PacBio
