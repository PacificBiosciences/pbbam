// File Description
/// \file ZmwReadStitcher.cpp
/// \brief Implements the ZmwReadStitcher class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/virtual/ZmwReadStitcher.h"

#include <deque>
#include <stdexcept>
#include <utility>

#include "VirtualZmwReader.h"
#include "pbbam/DataSet.h"
#include "pbbam/EntireFileQuery.h"
#include "pbbam/MakeUnique.h"
#include "pbbam/PbiFilter.h"
#include "pbbam/PbiFilterQuery.h"

namespace PacBio {
namespace BAM {

struct ZmwReadStitcher::ZmwReadStitcherPrivate
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
        // set up source queue
        std::string primaryFn;
        std::string scrapsFn;
        const auto& resources = dataset.ExternalResources();
        for (const ExternalResource& resource : resources) {

            primaryFn.clear();
            scrapsFn.clear();

            // if resource is possible "primary" BAM
            const auto& metatype = resource.MetaType();
            if (metatype == "PacBio.SubreadFile.SubreadBamFile" ||
                metatype == "PacBio.SubreadFile.HqRegionBamFile") {
                // possible resolve relative path
                primaryFn = dataset.ResolvePath(resource.ResourceId());

                // check for associated scraps file
                const auto& childResources = resource.ExternalResources();
                for (const auto& childResource : childResources) {
                    const auto& childMetatype = childResource.MetaType();
                    if (childMetatype == "PacBio.SubreadFile.ScrapsBamFile" ||
                        childMetatype == "PacBio.SubreadFile.HqScrapsBamFile") {
                        // possible resolve relative path
                        scrapsFn = dataset.ResolvePath(childResource.ResourceId());
                        break;
                    }
                }
            }

            // queue up source for later
            if (!primaryFn.empty() && !scrapsFn.empty())
                sources_.emplace_back(std::make_pair(primaryFn, scrapsFn));
        }

        OpenNextReader();
    }

public:
    bool HasNext() const { return (currentReader_ && currentReader_->HasNext()); }

    VirtualZmwBamRecord Next()
    {
        if (currentReader_) {
            const auto result = currentReader_->Next();
            if (!currentReader_->HasNext()) OpenNextReader();
            return result;
        }

        // no reader active
        throw std::runtime_error{
            "no readers active, make sure you use "
            "ZmwReadStitcher::HasNext before "
            "requesting next record"};
    }

    std::vector<BamRecord> NextRaw()
    {
        if (currentReader_) {
            const auto result = currentReader_->NextRaw();
            if (!currentReader_->HasNext()) OpenNextReader();
            return result;
        }

        // no reader active
        throw std::runtime_error{
            "no readers active, make sure you use "
            "ZmwReadStitcher::HasNext before "
            "requesting next group of records"};
    }

    BamHeader PrimaryHeader() const { return currentReader_->PrimaryHeader(); }

    BamHeader ScrapsHeader() const { return currentReader_->ScrapsHeader(); }

private:
    std::deque<std::pair<std::string, std::string> > sources_;
    std::unique_ptr<internal::VirtualZmwReader> currentReader_;
    PbiFilter filter_;

private:
    void OpenNextReader()
    {
        currentReader_.reset(nullptr);

        // find next source pair with data
        while (!sources_.empty()) {
            const auto nextSource = sources_.front();
            sources_.pop_front();

            currentReader_ = std::make_unique<internal::VirtualZmwReader>(
                nextSource.first, nextSource.second, filter_);
            if (currentReader_->HasNext()) return;
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

ZmwReadStitcher::~ZmwReadStitcher() {}

bool ZmwReadStitcher::HasNext() { return d_->HasNext(); }

VirtualZmwBamRecord ZmwReadStitcher::Next() { return d_->Next(); }

std::vector<BamRecord> ZmwReadStitcher::NextRaw() { return d_->NextRaw(); }

BamHeader ZmwReadStitcher::PrimaryHeader() const { return d_->PrimaryHeader().DeepCopy(); }

BamHeader ZmwReadStitcher::ScrapsHeader() const { return d_->ScrapsHeader().DeepCopy(); }

}  // namespace BAM
}  // namespace PacBio
