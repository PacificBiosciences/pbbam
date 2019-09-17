// File Description
/// \file BamFileMerger.cpp
/// \brief Implements the BamFileMerger & helper classes.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/BamFileMerger.h"

#include <memory>
#include <set>
#include <sstream>
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

using CompositeMergeItem = internal::CompositeMergeItem;

struct QNameSorter
{
    bool operator()(const CompositeMergeItem& lhs, const CompositeMergeItem& rhs) const
    {
        const BamRecord& l = lhs.record;
        const BamRecord& r = rhs.record;

        // movie name
        std::string lMovieName, rMovieName;  // TODO(CD): memoize movienames?
        try {
            lMovieName = l.MovieName();
        } catch (std::runtime_error const& err) {
            std::ostringstream s;
            s << "[pbbam] BAM file merging ERROR: could not get movie name:\n"
              << "  file: " << lhs.reader->Filename() << '\n'
              << "  reason: " << err.what();
            throw std::runtime_error(s.str());
        }
        try {
            rMovieName = r.MovieName();
        } catch (std::runtime_error const& err) {
            std::ostringstream s;
            s << "[pbbam] BAM file merging ERROR: could not get movie name:\n"
              << "  file: " << rhs.reader->Filename() << '\n'
              << "  reason: " << err.what();
            throw std::runtime_error(s.str());
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
        if (lhsQStart != rhsQStart) return lhsQStart < rhsQStart;

        const auto lhsQEnd = l.QueryEnd();
        const auto rhsQEnd = r.QueryEnd();
        return lhsQEnd < rhsQEnd;
    }
};

class ICollator
{
public:
    virtual bool GetNext(BamRecord&) = 0;
    virtual ~ICollator() = default;

protected:
    ICollator() = default;
};

template <typename Comp>
class CollatorImpl : public ICollator
{
public:
    CollatorImpl(std::vector<std::unique_ptr<BamReader>> readers) : ICollator()
    {
        for (auto&& reader : readers) {
            auto item = CompositeMergeItem{std::move(reader)};
            if (item.reader->GetNext(item.record)) mergeItems_.insert(std::move(item));
        }
    }

    bool GetNext(BamRecord& record) override
    {
        if (mergeItems_.empty()) return false;

        // Move first record into our result
        auto& firstItem = const_cast<CompositeMergeItem&>(*mergeItems_.begin());
        auto& firstRecord = firstItem.record;
        std::swap(record, firstRecord);

        // Try to read next record from current reader. If available, re-insert
        // into the set. Otherwise, just drop it (dtor will release resource).
        CompositeMergeItem tmp(std::move(firstItem));
        mergeItems_.erase(mergeItems_.begin());
        if (tmp.reader->GetNext(tmp.record)) mergeItems_.insert(std::move(tmp));
        return true;
    }

private:
    std::multiset<CompositeMergeItem, Comp> mergeItems_;
};

using QNameCollator = CollatorImpl<QNameSorter>;
using AlignedCollator = CollatorImpl<PositionSorter>;

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
    std::unique_ptr<ICollator> collator;
    if (isCoordinateSorted)
        collator = std::make_unique<AlignedCollator>(std::move(readers));
    else
        collator = std::make_unique<QNameCollator>(std::move(readers));
    return collator;
}

std::unique_ptr<IRecordWriter> MakeBamWriter(const std::vector<std::unique_ptr<BamReader>>& readers,
                                             const std::string& outputFilename,
                                             const bool createPbi, const ProgramInfo& pgInfo)
{
    if (outputFilename.empty())
        throw std::runtime_error{"[pbbam] BAM file merging ERROR: no output filename provided"};

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
            throw std::runtime_error{
                "[pbbam] BAM file merging ERROR: BAM file sort orders do not match, aborting "
                "merge"};
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

}  // namespace

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
    if (bamFiles.empty())
        throw std::runtime_error{"[pbbam] BAM file merging ERROR: no input filenames provided"};

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
    if (bamFiles.empty())
        throw std::runtime_error{"[pbbam] BAM file merging ERROR: no input filenames provided"};

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
    if (bamFiles.empty())
        throw std::runtime_error{"[pbbam] BAM file merging ERROR: no input filenames provided"};

    auto readers = MakeBamReaders(std::move(bamFiles), PbiFilter::FromDataSet(dataset));
    const bool isCoordinateSorted = readers.front()->Header().SortOrder() == "coordinate";

    auto collator = MakeCollator(std::move(readers), isCoordinateSorted);

    BamRecord record;
    while (collator->GetNext(record))
        writer.Write(record);
}

}  // namespace BAM
}  // namespace PacBio
