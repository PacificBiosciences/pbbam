// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/internal/DataSetBaseTypes.h"

#include <cstddef>

#include <boost/algorithm/string.hpp>

#include "DataSetUtils.h"
#include "TimeUtils.h"
#include "pbbam/DataSetTypes.h"

namespace PacBio {
namespace BAM {
namespace internal {

// ----------------
// BaseEntityType
// ----------------

BaseEntityType::BaseEntityType(const std::string& label, const XsdType& xsd)
    : DataSetElement(label, xsd)
{
    if (Version().empty()) Version(internal::XML_VERSION);
}

DEFINE_ACCESSORS(BaseEntityType, Extensions, Extensions)

BaseEntityType& BaseEntityType::Extensions(const PacBio::BAM::Extensions& extensions)
{
    Extensions() = extensions;
    return *this;
}

// ----------------
// DataEntityType
// ----------------

DataEntityType::DataEntityType(const std::string& label, const XsdType& xsd)
    : BaseEntityType(label, xsd)
{
}

// -----------------
// IndexedDataType
// -----------------

IndexedDataType::IndexedDataType(const std::string& metatype, const std::string& filename,
                                 const std::string& label, const XsdType& xsd)
    : InputOutputDataType(metatype, filename, label, xsd)
{
}

DEFINE_ACCESSORS(IndexedDataType, FileIndices, FileIndices)

IndexedDataType& IndexedDataType::FileIndices(const PacBio::BAM::FileIndices& indices)
{
    FileIndices() = indices;
    return *this;
}

// ---------------------
// InputOutputDataType
// ---------------------

InputOutputDataType::InputOutputDataType(const std::string& metatype, const std::string& filename,
                                         const std::string& label, const XsdType& xsd)
    : StrictEntityType(metatype, label, xsd)
{
    ResourceId(filename);
}

// ----------------
// StrictEntityType
// ----------------

StrictEntityType::StrictEntityType(const std::string& metatype, const std::string& label,
                                   const XsdType& xsd)
    : BaseEntityType(label, xsd)
{
    // MetaType
    MetaType(metatype);

    // TimeStampedName
    const size_t numChars = metatype.size();
    std::string transformedMetatype;
    transformedMetatype.resize(numChars);
    for (size_t i = 0; i < numChars; ++i) {
        const char c = metatype.at(i);
        transformedMetatype[i] = ((c == '.') ? '_' : tolower(c));
    }
    const std::string tsn =
        transformedMetatype + "-" + internal::ToDataSetFormat(internal::CurrentTime());
    TimeStampedName(tsn);

    // UniqueId
    UniqueId(internal::GenerateUuid());
}

}  // namespace internal
}  // namespace BAM
}  // namespace PacBio
