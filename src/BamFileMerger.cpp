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
#include "pbbam/IRecordWriter.h"
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
        std::string lMovieName, rMovieName; // TODO(CD): memoize movienames?
        try {
            lMovieName = l.MovieName();
        } catch (std::runtime_error const& err) {
            std::string msg{lhs.reader->Filename() + ": Could not get MovieName for BamFile. " + err.what()};
            throw std::runtime_error(msg);
        }
        try {
            rMovieName = r.MovieName();
        } catch (std::runtime_error const& err) {
            std::string msg{lhs.reader->Filename() + ": Could not get MovieName for BamFile. " + err.what()};
            throw std::runtime_error(msg);
        }
        const int cmp = lMovieName.compare(rMovieName);
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

std::vector<std::unique_ptr<BamReader>> MakeBamReaders(std::vector<BamFile> bamFiles,
                                                       PbiFilter filter = PbiFilter{})
{
    std::vector<std::unique_ptr<BamReader>> readers;
    for (auto& file : bamFiles) {
        if (filter.IsEmpty())
            readers.emplace_back(std::make_unique<BamReader>(std::move(file)));
        else
            readers.emplace_back(std::make_unique<PbiIndexedBamReader>(filter, std::move(file)));
    }
    assert(!readers.empty());
    return readers;
}

std::unique_ptr<ICollator> MakeCollator(std::vector<std::unique_ptr<BamReader>> readers,
                                        const bool isCoordinateSorted = false)
{
    if (isCoordinateSorted)
        return std::make_unique<AlignedCollator>(std::move(readers));
    else
        return std::make_unique<QNameCollator>(std::move(readers));
}

std::unique_ptr<IRecordWriter> MakeBamWriter(const std::vector<std::unique_ptr<BamReader>>& readers,
                                             const std::string& outputFilename,
                                             const bool createPbi, const ProgramInfo& pgInfo)
{
    if (outputFilename.empty())
        throw std::runtime_error{"no output BAM filename provide to BamFileMerger"};

    // read headers
    std::vector<BamHeader> headers;
    for (const auto& reader : readers)
        headers.push_back(reader->Header());
    assert(!headers.empty());

    // merge headers
    BamHeader mergedHeader = headers.at(0);
    const std::string usingSortOrder = mergedHeader.SortOrder();
    for (size_t i = 1; i < headers.size(); ++i) {
        const auto& header = headers.at(i);
        if (header.SortOrder() != usingSortOrder)
            throw std::runtime_error{"BAM file sort orders do not match, aborting merge"};
        mergedHeader += header;
    }

    // maybe add program info
    if (pgInfo.IsValid()) mergedHeader.AddProgram(pgInfo);

    // create BAM writer (PBI-on-the-fly?)
    if (createPbi)
        return std::make_unique<IndexedBamWriter>(outputFilename, mergedHeader);
    else
        return std::make_unique<BamWriter>(outputFilename, mergedHeader);
}

}  // namespace anonymous

void BamFileMerger::Merge(const std::vector<std::string>& bamFilenames,
                          const std::string& outputFilename, bool createPbi,
                          const ProgramInfo& pgInfo)
{
    std::vector<BamFile> bamFiles;
    for (const auto& fn : bamFilenames)
        bamFiles.emplace_back(fn);

    auto readers = MakeBamReaders(std::move(bamFiles));
    const bool isCoordinateSorted = readers.front()->Header().SortOrder() == "coordinate";

    auto writer = MakeBamWriter(readers, outputFilename, createPbi, pgInfo);
    auto collator = MakeCollator(std::move(readers), isCoordinateSorted);

    BamRecord record;
    while (collator->GetNext(record))
        writer->Write(record);
}

void BamFileMerger::Merge(const DataSet& dataset, const std::string& outputFilename, bool createPbi,
                          const ProgramInfo& pgInfo)
{
    std::vector<BamFile> bamFiles = dataset.BamFiles();
    if (bamFiles.empty()) throw std::runtime_error{"no input filenames provided to BamFileMerger"};

    auto readers = MakeBamReaders(std::move(bamFiles), PbiFilter::FromDataSet(dataset));
    const bool isCoordinateSorted = readers.front()->Header().SortOrder() == "coordinate";

    auto writer = MakeBamWriter(readers, outputFilename, createPbi, pgInfo);
    auto collator = MakeCollator(std::move(readers), isCoordinateSorted);

    BamRecord record;
    while (collator->GetNext(record))
        writer->Write(record);
}

void BamFileMerger::Merge(const std::vector<std::string>& bamFilenames, IRecordWriter& writer)
{
    std::vector<BamFile> bamFiles;
    for (const auto& fn : bamFilenames)
        bamFiles.emplace_back(fn);
    if (bamFiles.empty()) throw std::runtime_error{"no input filenames provided to BamFileMerger"};

    auto readers = MakeBamReaders(std::move(bamFiles));
    const bool isCoordinateSorted = readers.front()->Header().SortOrder() == "coordinate";

    auto collator = MakeCollator(std::move(readers), isCoordinateSorted);

    BamRecord record;
    while (collator->GetNext(record))
        writer.Write(record);
}

void BamFileMerger::Merge(const DataSet& dataset, IRecordWriter& writer)
{
    std::vector<BamFile> bamFiles = dataset.BamFiles();
    if (bamFiles.empty()) throw std::runtime_error{"no input filenames provided to BamFileMerger"};

    auto readers = MakeBamReaders(std::move(bamFiles), PbiFilter::FromDataSet(dataset));
    const bool isCoordinateSorted = readers.front()->Header().SortOrder() == "coordinate";

    auto collator = MakeCollator(std::move(readers), isCoordinateSorted);

    BamRecord record;
    while (collator->GetNext(record))
        writer.Write(record);
}

}  // namespace BAM
}  // namespace PacBio
