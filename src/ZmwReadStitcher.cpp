#include "PbbamInternalConfig.h"

#include <pbbam/virtual/ZmwReadStitcher.h>

#include <deque>
#include <stdexcept>
#include <utility>

#include <pbbam/DataSet.h>
#include <pbbam/EntireFileQuery.h>
#include <pbbam/PbiFilter.h>
#include <pbbam/PbiFilterQuery.h>

#include "VirtualStitching.h"
#include "VirtualZmwReader.h"

namespace PacBio {
namespace BAM {

class ZmwReadStitcher::ZmwReadStitcherPrivate
{
public:
    ZmwReadStitcherPrivate(std::string primaryBamFilePath, std::string scrapsBamFilePath,
                           PbiFilter filter)
        : filter_{std::move(filter)}
    {
        sources_.push_back({std::move(primaryBamFilePath), std::move(scrapsBamFilePath)});
        OpenNextReader();
    }

    ZmwReadStitcherPrivate(const DataSet& dataset) : filter_{PbiFilter::FromDataSet(dataset)}
    {
        sources_ = SourcesFromDataset(dataset);
        OpenNextReader();
    }

    bool HasNext() const { return (currentReader_ && currentReader_->HasNext()); }

    VirtualZmwBamRecord Next()
    {
        if (currentReader_) {
            const auto result = currentReader_->Next();
            if (!currentReader_->HasNext()) {
                OpenNextReader();
            }
            return result;
        }

        // no reader active
        throw std::runtime_error{
            "[pbbam] ZMW stitching ERROR: "
            "no readers active, make sure you use "
            "ZmwReadStitcher::HasNext before "
            "requesting next record"};
    }

    std::vector<BamRecord> NextRaw()
    {
        if (currentReader_) {
            const auto result = currentReader_->NextRaw();
            if (!currentReader_->HasNext()) {
                OpenNextReader();
            }
            return result;
        }

        // no reader active
        throw std::runtime_error{
            "[pbbam] ZMW stitching ERROR: "
            "no readers active, make sure you use "
            "ZmwReadStitcher::HasNext before "
            "requesting next group of records"};
    }

    BamHeader PrimaryHeader() const { return currentReader_->PrimaryHeader(); }

    BamHeader ScrapsHeader() const { return currentReader_->ScrapsHeader(); }

    BamHeader StitchedHeader() const { return currentReader_->StitchedHeader(); }

private:
    StitchingSources sources_;
    std::unique_ptr<VirtualZmwReader> currentReader_;
    PbiFilter filter_;

    void OpenNextReader()
    {
        currentReader_.reset(nullptr);

        // find next source pair with data
        while (!sources_.empty()) {
            const auto nextSource = sources_.front();
            sources_.pop_front();

            currentReader_ =
                std::make_unique<VirtualZmwReader>(nextSource.first, nextSource.second, filter_);
            if (currentReader_->HasNext()) {
                return;
            }
        }
    }
};

// --------------------------------
// ZmwReadStitcher implementation
// --------------------------------

ZmwReadStitcher::ZmwReadStitcher(std::string primaryBamFilePath, std::string scrapsBamFilePath)
    : ZmwReadStitcher{std::move(primaryBamFilePath), std::move(scrapsBamFilePath), PbiFilter{}}
{
}

ZmwReadStitcher::ZmwReadStitcher(std::string primaryBamFilePath, std::string scrapsBamFilePath,
                                 PbiFilter filter)
    : d_{std::make_unique<ZmwReadStitcherPrivate>(std::move(primaryBamFilePath),
                                                  std::move(scrapsBamFilePath), std::move(filter))}
{
}

ZmwReadStitcher::ZmwReadStitcher(const DataSet& dataset)
    : d_{std::make_unique<ZmwReadStitcherPrivate>(dataset)}
{
}

ZmwReadStitcher::ZmwReadStitcher(ZmwReadStitcher&&) noexcept = default;

ZmwReadStitcher& ZmwReadStitcher::operator=(ZmwReadStitcher&&) noexcept = default;

ZmwReadStitcher::~ZmwReadStitcher() = default;

bool ZmwReadStitcher::HasNext() { return d_->HasNext(); }

VirtualZmwBamRecord ZmwReadStitcher::Next() { return d_->Next(); }

std::vector<BamRecord> ZmwReadStitcher::NextRaw() { return d_->NextRaw(); }

BamHeader ZmwReadStitcher::PrimaryHeader() const { return d_->PrimaryHeader().DeepCopy(); }

BamHeader ZmwReadStitcher::ScrapsHeader() const { return d_->ScrapsHeader().DeepCopy(); }

BamHeader ZmwReadStitcher::StitchedHeader() const { return d_->StitchedHeader().DeepCopy(); }

}  // namespace BAM
}  // namespace PacBio
