#ifndef PBBAM_COMPOSITEREADER_INL
#define PBBAM_COMPOSITEREADER_INL

#include <pbbam/Config.h>

#include <pbbam/CompositeBamReader.h>

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
bool CompositeMergeItemSorter<CompareType>::operator()(const CompositeMergeItem& lhs,
                                                       const CompositeMergeItem& rhs) const
{
    const auto& l = lhs.record;
    const auto& r = rhs.record;
    return CompareType()(l, r);
}

}  // namespace internal

// -----------------------------------
// general SortedCompositeBamReader
// -----------------------------------

template <typename OrderByType>
SortedCompositeBamReader<OrderByType>::SortedCompositeBamReader(const DataSet& dataset)
    : SortedCompositeBamReader(dataset.BamFiles())
{
}

template <typename OrderByType>
SortedCompositeBamReader<OrderByType>::SortedCompositeBamReader(std::vector<BamFile> bamFiles)
    : internal::IQuery{}, bamFiles_{std::move(bamFiles)}
{
    // create readers for files
    for (const auto& bamFile : bamFiles_) {
        internal::CompositeMergeItem item{std::make_unique<BamReader>(bamFile)};
        if (item.reader->GetNext(item.record)) {
            mergeItems_.insert(std::move(item));
        }
    }
}

template <typename OrderByType>
SortedCompositeBamReader<OrderByType>::SortedCompositeBamReader(
    SortedCompositeBamReader&&) noexcept = default;

template <typename OrderByType>
SortedCompositeBamReader<OrderByType>& SortedCompositeBamReader<OrderByType>::operator=(
    SortedCompositeBamReader&&) noexcept = default;

template <typename OrderByType>
SortedCompositeBamReader<OrderByType>::~SortedCompositeBamReader() = default;

template <typename OrderByType>
bool SortedCompositeBamReader<OrderByType>::GetNext(BamRecord& record)
{
    if (mergeItems_.empty()) {
        return false;
    }

    // Move first record into our result
    auto& firstItem = const_cast<internal::CompositeMergeItem&>(*mergeItems_.begin());
    auto& firstRecord = firstItem.record;
    std::swap(record, firstRecord);

    // Try to read next record from current reader. If available, re-insert
    // into the set. Otherwise, just drop it (dtor will release resource).
    internal::CompositeMergeItem tmp(std::move(firstItem));
    mergeItems_.erase(mergeItems_.begin());
    if (tmp.reader->GetNext(tmp.record)) {
        mergeItems_.insert(std::move(tmp));
    }
    return true;
}

// ------------------------------
// PbiFilterCompositeReader
// ------------------------------

template <typename OrderByType>
PbiFilterCompositeBamReader<OrderByType>::PbiFilterCompositeBamReader(
    const PbiFilter& filter, const std::vector<BamFile>& bamFiles)
    : PbiFilterCompositeBamReader<OrderByType>{filter, bamFiles, MakePbiIndexCache(bamFiles)}
{
}

template <typename OrderByType>
PbiFilterCompositeBamReader<OrderByType>::PbiFilterCompositeBamReader(
    const PbiFilter& filter, const std::vector<BamFile>& bamFiles, const PbiIndexCache& cache)
    : SortedCompositeBamReader<OrderByType>{bamFiles}, indexCache_{cache}, numReads_{0}
{
    Filter(filter);
}

template <typename OrderByType>
PbiFilterCompositeBamReader<OrderByType>::PbiFilterCompositeBamReader(const PbiFilter& filter,
                                                                      const DataSet& dataset)
    : PbiFilterCompositeBamReader<OrderByType>{filter, dataset.BamFiles()}
{
}

template <typename OrderByType>
PbiFilterCompositeBamReader<OrderByType>::PbiFilterCompositeBamReader(const PbiFilter& filter,
                                                                      const DataSet& dataset,
                                                                      const PbiIndexCache& cache)
    : PbiFilterCompositeBamReader<OrderByType>{filter, dataset.BamFiles(), cache}
{
}

template <typename OrderByType>
PbiFilterCompositeBamReader<OrderByType>& PbiFilterCompositeBamReader<OrderByType>::Filter(
    const PbiFilter& filter)
{
    // reset reader queue
    this->mergeItems_.clear();

    // create readers for files
    container_type updatedMergeItems;
    std::vector<std::string> missingPbi;
    for (size_t i = 0; i < this->bamFiles_.size(); ++i) {
        const auto& bamFile = this->bamFiles_.at(i);
        if (bamFile.PacBioIndexExists()) {
            auto item = internal::CompositeMergeItem{std::unique_ptr<BamReader>{
                new PbiIndexedBamReader{filter, std::move(bamFile), indexCache_->at(i)}}};
            if (item.reader->GetNext(item.record)) {
                updatedMergeItems.insert(std::move(item));
            }
            // else not an error, simply no data matching filter
        } else {
            missingPbi.push_back(bamFile.Filename());
        }
    }

    // throw if any files missing PBI
    if (!missingPbi.empty()) {
        std::ostringstream e;
        e << "[pbbam] composite BAM reader ERROR: failed to open because the following files are "
             "missing a *.pbi index:\n";
        for (const auto& fn : missingPbi) {
            e << "  " << fn << '\n';
        }
        throw std::runtime_error{e.str()};
    }

    // update our actual container, store num matching reads, sort & and return
    this->mergeItems_ = std::move(updatedMergeItems);
    numReads_ = 0;
    for (const auto& item : this->mergeItems_) {
        auto* pbiReader = dynamic_cast<PbiIndexedBamReader*>(item.reader.get());
        numReads_ += pbiReader->NumReads();
    }
    return *this;
}

template <typename OrderByType>
uint32_t PbiFilterCompositeBamReader<OrderByType>::NumReads() const
{
    return numReads_;
}

}  // namespace BAM
}  // namespace PacBio

#endif // PBBAM_COMPOSITEREADER_INL
