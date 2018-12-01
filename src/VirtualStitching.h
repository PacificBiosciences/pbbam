// File Description
/// \file VirtualStitching.h
/// \brief Utilities for virtual ZMW stitiching.
//
// Author: Derek Barnett

#ifndef VIRTUALSTITCHING_H
#define VIRTUALSTITCHING_H

#include <deque>
#include <string>
#include <utility>

#include <boost/optional.hpp>

#include "pbbam/DataSet.h"

namespace PacBio {
namespace BAM {

using StitchingSources = std::deque<std::pair<std::string, std::string>>;

inline boost::optional<std::string> ScrapsFileId(const ExternalResource& resource)
{
    const auto& childResources = resource.ExternalResources();
    for (const auto& childResource : childResources) {
        const auto& childMetatype = childResource.MetaType();
        if (childMetatype == "PacBio.SubreadFile.ScrapsBamFile" ||
            childMetatype == "PacBio.SubreadFile.HqScrapsBamFile") {
            return childResource.ResourceId();
        }
    }
    return boost::none;
}

inline StitchingSources SourcesFromDataset(const DataSet& dataset)
{
    StitchingSources sources;

    const ExternalResources& resources = dataset.ExternalResources();
    for (const ExternalResource& resource : resources) {

        boost::optional<std::string> primaryId;
        boost::optional<std::string> scrapsId;

        // if resource is possible "primary" BAM, store & look for associated scraps
        const auto& metatype = resource.MetaType();
        if (metatype == "PacBio.SubreadFile.SubreadBamFile" ||
            metatype == "PacBio.SubreadFile.HqRegionBamFile") {
            primaryId = resource.ResourceId();
            scrapsId = ScrapsFileId(resource);
        }

        // if found, resolve paths & store
        if (primaryId && scrapsId) {
            std::string primaryFn = dataset.ResolvePath(primaryId.get());
            std::string scrapsFn = dataset.ResolvePath(scrapsId.get());
            sources.emplace_back(std::make_pair(primaryFn, scrapsFn));
        }
    }

    return sources;
}

}  // namespace BAM
}  // namespace PacBio

#endif  // VIRTUALSTITCHING_H