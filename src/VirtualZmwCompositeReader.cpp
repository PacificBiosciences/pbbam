// File Description
/// \file VirtualZmwCompositeReader.cpp
/// \brief Implements the VirtualZmwCompositeReader class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "VirtualZmwCompositeReader.h"

#include <boost/algorithm/string.hpp>

#include "pbbam/MakeUnique.h"

namespace PacBio {
namespace BAM {
namespace internal {

VirtualZmwCompositeReader::VirtualZmwCompositeReader(const DataSet& dataset)
    : currentReader_(nullptr), filter_(PbiFilter::FromDataSet(dataset))
{
    // set up source queue
    std::string primaryFn;
    std::string scrapsFn;
    const ExternalResources& resources = dataset.ExternalResources();
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

    // open first available source
    OpenNextReader();
}

bool VirtualZmwCompositeReader::HasNext() { return (currentReader_ && currentReader_->HasNext()); }

VirtualZmwBamRecord VirtualZmwCompositeReader::Next()
{
    if (currentReader_) {
        const auto result = currentReader_->Next();
        if (!currentReader_->HasNext()) OpenNextReader();
        return result;
    }

    // no reader active
    throw std::runtime_error{
        "no readers active, make sure you use "
        "VirtualZmwCompositeReader::HasNext before "
        "requesting next record"};
}

std::vector<BamRecord> VirtualZmwCompositeReader::NextRaw()
{
    if (currentReader_) {
        const auto result = currentReader_->NextRaw();
        if (!currentReader_->HasNext()) OpenNextReader();
        return result;
    }

    // no reader active
    throw std::runtime_error{
        "no readers active, make sure you use "
        "VirtualZmwCompositeReader::HasNext before "
        "requesting next group of records"};
}

void VirtualZmwCompositeReader::OpenNextReader()
{
    currentReader_.reset(nullptr);

    // find next source pair with data
    while (!sources_.empty()) {
        const auto nextSource = sources_.front();
        sources_.pop_front();

        currentReader_ =
            std::make_unique<VirtualZmwReader>(nextSource.first, nextSource.second, filter_);
        if (currentReader_->HasNext()) return;
    }
}

}  // namespace internal
}  // namespace BAM
}  // namespace PacBio
