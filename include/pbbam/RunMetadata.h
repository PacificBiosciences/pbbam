#ifndef PBBAM_RUNMETADATA_H
#define PBBAM_RUNMETADATA_H

#include <pbbam/Config.h>

#include <pbbam/CollectionMetadata.h>

#include <iosfwd>
#include <map>
#include <string>

namespace PacBio {
namespace BAM {

class RunMetadata
{
public:
    ///
    /// \brief Read a single CollectionMetadata from XML file (or input stream).
    ///
    /// Intended for use with '<id>.metadata.xml'.
    ///
    /// \throw std::runtime_error on failure to parse XML for CollectionMetadata,
    ///                           or if the number of SubreadSets found is != 1.
    ///
    static CollectionMetadata Collection(const std::string& metadataXmlFn);
    static CollectionMetadata Collection(std::istream& in);

    ///
    /// \brief Read multiple CollectionMetadata objects from XML file (or input stream).
    ///
    /// Intended for use with multi-collection 'id.run.metadata.xml', but will
    /// work for single collection input.
    ///
    /// \throw std::runtime_error on failure to parse XML for CollectionMetadata(s)
    ///
    /// \returns map {subread set name => CollectionMetadata}
    ///
    static std::map<std::string, CollectionMetadata> Collections(
        const std::string& runMetadataXmlFn);
    static std::map<std::string, CollectionMetadata> Collections(std::istream& in);
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_RUNMETADATA_H
