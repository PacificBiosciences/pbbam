// File Description
/// \file BamFileMerger.cpp
/// \brief Implements the BamFileMerger & helper classes.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/BamFileMerger.h"

#include <deque>
#include <memory>
#include <stdexcept>
#include <vector>

#include "pbbam/BamFile.h"
#include "pbbam/BamHeader.h"
#include "pbbam/BamReader.h"
#include "pbbam/BamRecord.h"
#include "pbbam/CompositeBamReader.h"
#include "pbbam/DataSet.h"
#include "pbbam/IndexedBamWriter.h"
#include "pbbam/PbiBuilder.h"
#include "pbbam/PbiFilter.h"
#include "pbbam/PbiIndexedBamReader.h"
#include "pbbam/RecordType.h"

namespace PacBio {
namespace BAM {
namespace {  // anonymous

class ICollator
{
public:
    virtual ~ICollator(void) = default;

    bool GetNext(BamRecord& record)
    {
        // nothing left to read
        if (mergeItems_.empty()) return false;

        // non-destructive 'pop' of first item from queue
        auto firstIter = mergeItems_.begin();
        auto firstItem = internal::CompositeMergeItem{std::move(firstIter->reader),
                                                      std::move(firstIter->record)};
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

protected:
    std::deque<internal::CompositeMergeItem> mergeItems_;

protected:
    ICollator(std::vector<std::unique_ptr<PacBio::BAM::BamReader>> readers)
    {
        for (auto&& reader : readers) {
            auto item = internal::CompositeMergeItem{std::move(reader)};
            if (item.reader->GetNext(item.record)) mergeItems_.push_back(std::move(item));
        }
    }

    virtual void UpdateSort(void) = 0;
};

struct QNameSorter
    : std::binary_function<internal::CompositeMergeItem, internal::CompositeMergeItem, bool>
{
    bool operator()(const internal::CompositeMergeItem& lhs,
                    const internal::CompositeMergeItem& rhs)
    {
        const BamRecord& l = lhs.record;
        const BamRecord& r = rhs.record;

        // movie name
        const int cmp = l.MovieName().compare(r.MovieName());
        if (cmp != 0) return cmp < 0;

        // hole number
        const auto lhsZmw = l.HoleNumber();
        const auto rhsZmw = r.HoleNumber();
        if (lhsZmw != rhsZmw) return lhsZmw < rhsZmw;

        // shuffle CCS/transcript reads after all others
        if (IsCcsOrTranscript(l.Type())) return false;
        if (IsCcsOrTranscript(r.Type())) return true;

        // sort on qStart, then finally qEnd
        const auto lhsQStart = l.QueryStart();
        const auto rhsQStart = r.QueryStart();
        return lhsQStart < rhsQStart;
    }
};

class QNameCollator : public ICollator
{
public:
    QNameCollator(std::vector<std::unique_ptr<PacBio::BAM::BamReader>> readers)
        : ICollator(std::move(readers))
    {
        UpdateSort();
    }

    void UpdateSort(void) { std::sort(mergeItems_.begin(), mergeItems_.end(), QNameSorter{}); }
};

class AlignedCollator : public ICollator
{
public:
    AlignedCollator(std::vector<std::unique_ptr<BamReader>> readers) : ICollator(std::move(readers))
    {
        UpdateSort();
    }

    void UpdateSort(void) { std::sort(mergeItems_.begin(), mergeItems_.end(), PositionSorter{}); }
};

void MergeImpl(std::vector<BamFile> bamFiles, const std::string& outputFilename,
               const PbiFilter& filter, bool createPbi, BamHeader initialOutputHeader)
{
    // I/O filenames check
    if (bamFiles.empty()) throw std::runtime_error{"no input filenames provided to BamFileMerger"};
    if (outputFilename.empty())
        throw std::runtime_error{"no output filename provide to BamFileMerger"};

    // attempt open input files
    std::vector<std::unique_ptr<BamReader>> readers;
    for (const auto& file : bamFiles) {
        if (filter.IsEmpty())
            readers.emplace_back(std::make_unique<BamReader>(file));
        else
            readers.emplace_back(std::make_unique<PbiIndexedBamReader>(filter, file));
    }
    assert(!readers.empty());

    // read headers
    std::vector<BamHeader> headers;
    for (const auto& reader : readers)
        headers.push_back(reader->Header());
    assert(!headers.empty());

    // merge headers
    BamHeader mergedHeader = initialOutputHeader;
    const std::string usingSortOrder = mergedHeader.SortOrder();
    const bool isCoordinateSorted = (usingSortOrder == "coordinate");
    for (const auto& header : headers) {
        if (header.SortOrder() != usingSortOrder)
            throw std::runtime_error{"BAM file sort orders do not match, aborting merge"};
        mergedHeader += header;
    }

    // setup collator - sort order?
    //
    // NOTE: collator takes ownership of readers
    //
    std::unique_ptr<ICollator> collator;
    if (isCoordinateSorted)
        collator = std::make_unique<AlignedCollator>(std::move(readers));
    else
        collator = std::make_unique<QNameCollator>(std::move(readers));

    // setup writer - PBI-on-the-fly?
    std::unique_ptr<IRecordWriter> writer;
    if (createPbi)
        writer = std::make_unique<IndexedBamWriter>(outputFilename, mergedHeader);
    else
        writer = std::make_unique<BamWriter>(outputFilename, mergedHeader);

    // write collated records
    BamRecord record;
    while (collator->GetNext(record))
        writer->Write(record);
}

}  // namespace anonymous

void BamFileMerger::Merge(const std::vector<std::string>& bamFilenames,
                          const std::string& outputFilename, bool createPbi,
                          BamHeader initialOutputHeader)
{
    std::vector<BamFile> bamFiles;
    for (const auto& fn : bamFilenames)
        bamFiles.emplace_back(fn);

    MergeImpl(std::move(bamFiles), outputFilename, PbiFilter{}, createPbi, initialOutputHeader);
}

void BamFileMerger::Merge(const DataSet& dataset, const std::string& outputFilename, bool createPbi,
                          BamHeader initialOutputHeader)
{
    MergeImpl(dataset.BamFiles(), outputFilename, PbiFilter::FromDataSet(dataset), createPbi,
              initialOutputHeader);
}

}  // namespace BAM
}  // namespace PacBio
