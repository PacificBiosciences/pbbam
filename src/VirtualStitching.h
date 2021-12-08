#ifndef PBBAM_VIRTUALSTITCHING_H
#define PBBAM_VIRTUALSTITCHING_H

#include <pbbam/Config.h>

#include <pbbam/DataSet.h>

#include <boost/optional.hpp>

#include <deque>
#include <string>
#include <utility>

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
    return {};
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
            std::string primaryFn = dataset.ResolvePath(*primaryId);
            std::string scrapsFn = dataset.ResolvePath(*scrapsId);
            sources.emplace_back(std::make_pair(primaryFn, scrapsFn));
        }
    }

    return sources;
}

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_VIRTUALSTITCHING_H
