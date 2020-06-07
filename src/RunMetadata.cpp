// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "RunMetadataParser.h"
#include "pbbam/RunMetadata.h"

#include <iosfwd>

namespace PacBio {
namespace BAM {

// ----------------------
// RunMetadata
// ----------------------

CollectionMetadata RunMetadata::Collection(const std::string& metadataXmlFn)
{
    return RunMetadataParser::LoadCollection(metadataXmlFn);
}

CollectionMetadata RunMetadata::Collection(std::istream& in)
{
    return RunMetadataParser::LoadCollection(in);
}

std::map<std::string, CollectionMetadata> RunMetadata::Collections(
    const std::string& runMetadataXmlFn)
{
    return RunMetadataParser::LoadCollections(runMetadataXmlFn);
}
std::map<std::string, CollectionMetadata> RunMetadata::Collections(std::istream& in)
{
    return RunMetadataParser::LoadCollections(in);
}

}  // namespace BAM
}  // namespace PacBio
