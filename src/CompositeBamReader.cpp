#include "PbbamInternalConfig.h"

#include <pbbam/CompositeBamReader.h>

#include <sstream>
#include <stdexcept>

#include <pbbam/BamFile.h>

namespace PacBio {
namespace BAM {

// -----------------------------------
// GenomicIntervalCompositeBamReader
// -----------------------------------

GenomicIntervalCompositeBamReader::GenomicIntervalCompositeBamReader(
    const std::vector<BamFile>& bamFiles)
    : GenomicIntervalCompositeBamReader{bamFiles, MakeBaiIndexCache(bamFiles)}
{
}

GenomicIntervalCompositeBamReader::GenomicIntervalCompositeBamReader(
    const std::vector<BamFile>& bamFiles, const BaiIndexCache& cache)
    : SortedCompositeBamReader<Compare::AlignmentPosition>{bamFiles}, indexCache_{cache}
{
    // no interval set, so no valid readers
    mergeItems_.clear();
}

GenomicIntervalCompositeBamReader::GenomicIntervalCompositeBamReader(const DataSet& dataset)
    : GenomicIntervalCompositeBamReader{dataset.BamFiles()}
{
}

GenomicIntervalCompositeBamReader::GenomicIntervalCompositeBamReader(const DataSet& dataset,
                                                                     const BaiIndexCache& cache)
    : GenomicIntervalCompositeBamReader{dataset.BamFiles(), cache}
{
}

GenomicIntervalCompositeBamReader::GenomicIntervalCompositeBamReader(
    const GenomicInterval& interval, const std::vector<BamFile>& bamFiles)
    : GenomicIntervalCompositeBamReader{interval, bamFiles, MakeBaiIndexCache(bamFiles)}
{
}

GenomicIntervalCompositeBamReader::GenomicIntervalCompositeBamReader(
    const GenomicInterval& interval, const std::vector<BamFile>& bamFiles,
    const BaiIndexCache& cache)
    : GenomicIntervalCompositeBamReader{bamFiles, cache}
{
    Interval(interval);
}

GenomicIntervalCompositeBamReader::GenomicIntervalCompositeBamReader(
    const GenomicInterval& interval, const DataSet& dataset)
    : GenomicIntervalCompositeBamReader{interval, dataset.BamFiles()}
{
}

GenomicIntervalCompositeBamReader::GenomicIntervalCompositeBamReader(
    const GenomicInterval& interval, const DataSet& dataset, const BaiIndexCache& cache)
    : GenomicIntervalCompositeBamReader{interval, dataset.BamFiles(), cache}
{
}

const GenomicInterval& GenomicIntervalCompositeBamReader::Interval() const { return interval_; }

GenomicIntervalCompositeBamReader& GenomicIntervalCompositeBamReader::Interval(
    const GenomicInterval& interval)
{
    // reset readers
    mergeItems_.clear();

    // create readers for files
    container_type updatedMergeItems;
    std::vector<std::string> missingBai;
    for (size_t i = 0; i < bamFiles_.size(); ++i) {
        const auto& bamFile = bamFiles_.at(i);
        if (bamFile.StandardIndexExists()) {
            internal::CompositeMergeItem item{std::unique_ptr<BamReader>{
                new BaiIndexedBamReader{interval, bamFile, indexCache_->at(i)}}};
            if (item.reader->GetNext(item.record)) {
                updatedMergeItems.insert(std::move(item));
            }
            // else not an error, simply no data matching interval
        } else {
            // maybe handle PBI-backed interval searches if BAI missing, but for now treat as error
            missingBai.push_back(bamFile.Filename());
        }
    }

    // throw if any files missing BAI
    if (!missingBai.empty()) {
        std::ostringstream e;
        e << "[pbbam] composite BAM reader ERROR: failed to open because the following files are "
             "missing a *.bai index:\n";
        for (const auto& fn : missingBai) {
            e << "  " << fn << '\n';
        }
        throw std::runtime_error{e.str()};
    }

    // update our actual container and return
    mergeItems_ = std::move(updatedMergeItems);
    return *this;
}

// ------------------------------
// SequentialCompositeBamReader
// ------------------------------

SequentialCompositeBamReader::SequentialCompositeBamReader(std::vector<BamFile> bamFiles)
    : internal::IQuery{}
{
    for (const auto& bamFile : bamFiles) {
        readers_.emplace_back(std::make_unique<BamReader>(bamFile));
    }
}

SequentialCompositeBamReader::SequentialCompositeBamReader(const DataSet& dataset)
    : SequentialCompositeBamReader{dataset.BamFiles()}
{
}

bool SequentialCompositeBamReader::GetNext(BamRecord& record)
{
    // try first reader, if successful return true
    // else pop reader and try next, until all readers exhausted
    while (!readers_.empty()) {
        auto& reader = readers_.front();
        if (reader->GetNext(record)) {
            return true;
        } else {
            readers_.pop_front();
        }
    }

    // no readers available
    return false;
}

}  // namespace BAM
}  // namespace PacBio
