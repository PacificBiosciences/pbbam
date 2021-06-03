#include "PbbamInternalConfig.h"

#include <pbbam/internal/DataSetBaseTypes.h>

#include <cstddef>

#include <boost/algorithm/string.hpp>

#include <pbbam/DataSetTypes.h>

#include "DataSetUtils.h"
#include "TimeUtils.h"

namespace PacBio {
namespace BAM {
namespace internal {

// ----------------
// BaseEntityType
// ----------------

BaseEntityType::BaseEntityType(const std::string& label, const XsdType& xsd)
    : DataSetElement(label, xsd)
{
    if (CreatedAt().empty()) {
        CreatedAt(TimeUtils::ToIso8601(TimeUtils::CurrentTime()));
    }
    if (Version().empty()) {
        Version(XML_VERSION);
    }
}

BaseEntityType::BaseEntityType(const std::string& label, const FromInputXml& fromInputXml,
                               const XsdType& xsd)
    : DataSetElement(label, fromInputXml, xsd)
{
}

const std::string& BaseEntityType::CreatedAt() const { return Attribute("CreatedAt"); }

std::string& BaseEntityType::CreatedAt() { return Attribute("CreatedAt"); }

BaseEntityType& BaseEntityType::CreatedAt(const std::string& createdAt)
{
    Attribute("CreatedAt", createdAt);
    return *this;
}

const std::string& BaseEntityType::Description() const { return Attribute("Description"); }

std::string& BaseEntityType::Description() { return Attribute("Description"); }

BaseEntityType& BaseEntityType::Description(const std::string& description)
{
    Attribute("Description", description);
    return *this;
}

DEFINE_ACCESSORS(BaseEntityType, Extensions, Extensions)

BaseEntityType& BaseEntityType::Extensions(const PacBio::BAM::Extensions& extensions)
{
    Extensions() = extensions;
    return *this;
}

const std::string& BaseEntityType::Format() const { return Attribute("Format"); }

std::string& BaseEntityType::Format() { return Attribute("Format"); }

BaseEntityType& BaseEntityType::Format(const std::string& format)
{
    Attribute("Format", format);
    return *this;
}

const std::string& BaseEntityType::ModifiedAt() const { return Attribute("ModifiedAt"); }

std::string& BaseEntityType::ModifiedAt() { return Attribute("ModifiedAt"); }

BaseEntityType& BaseEntityType::ModifiedAt(const std::string& modifiedAt)
{
    Attribute("ModifiedAt", modifiedAt);
    return *this;
}

const std::string& BaseEntityType::Name() const { return Attribute("Name"); }

std::string& BaseEntityType::Name() { return Attribute("Name"); }

BaseEntityType& BaseEntityType::Name(const std::string& name)
{
    Attribute("Name", name);
    return *this;
}

const std::string& BaseEntityType::ResourceId() const { return Attribute("ResourceId"); }

std::string& BaseEntityType::ResourceId() { return Attribute("ResourceId"); }

BaseEntityType& BaseEntityType::ResourceId(const std::string& resourceId)
{
    Attribute("ResourceId", resourceId);
    return *this;
}

const std::string& BaseEntityType::Tags() const { return Attribute("Tags"); }

std::string& BaseEntityType::Tags() { return Attribute("Tags"); }

BaseEntityType& BaseEntityType::Tags(const std::string& tags)
{
    Attribute("Tags", tags);
    return *this;
}

const std::string& BaseEntityType::Version() const { return Attribute("Version"); }

std::string& BaseEntityType::Version() { return Attribute("Version"); }

BaseEntityType& BaseEntityType::Version(const std::string& version)
{
    Attribute("Version", version);
    return *this;
}

// ----------------
// DataEntityType
// ----------------

DataEntityType::DataEntityType(const std::string& label, const XsdType& xsd)
    : BaseEntityType(label, xsd)
{
}

DataEntityType::DataEntityType(const std::string& label, const FromInputXml& fromInputXml,
                               const XsdType& xsd)
    : BaseEntityType(label, fromInputXml, xsd)
{
}

const std::string& DataEntityType::Checksum() const { return ChildText("Checksum"); }

std::string& DataEntityType::Checksum() { return ChildText("Checksum"); }

DataEntityType& DataEntityType::Checksum(const std::string& checksum)
{
    ChildText("Checksum", checksum);
    return *this;
}

const std::string& DataEntityType::EncodedValue() const { return ChildText("EncodedValue"); }

std::string& DataEntityType::EncodedValue() { return ChildText("EncodedValue"); }

DataEntityType& DataEntityType::EncodedValue(const std::string& encodedValue)
{
    ChildText("EncodedValue", encodedValue);
    return *this;
}

const std::string& DataEntityType::MetaType() const { return Attribute("MetaType"); }

std::string& DataEntityType::MetaType() { return Attribute("MetaType"); }

DataEntityType& DataEntityType::MetaType(const std::string& metatype)
{
    Attribute("MetaType", metatype);
    return *this;
}

const std::string& DataEntityType::SimpleValue() const { return Attribute("SimpleValue"); }

std::string& DataEntityType::SimpleValue() { return Attribute("SimpleValue"); }

DataEntityType& DataEntityType::SimpleValue(const std::string& simpleValue)
{
    Attribute("SimpleValue", simpleValue);
    return *this;
}

const std::string& DataEntityType::TimeStampedName() const { return Attribute("TimeStampedName"); }

std::string& DataEntityType::TimeStampedName() { return Attribute("TimeStampedName"); }

DataEntityType& DataEntityType::TimeStampedName(const std::string& timeStampedName)
{
    Attribute("TimeStampedName", timeStampedName);
    return *this;
}

const std::string& DataEntityType::UniqueId() const { return Attribute("UniqueId"); }

std::string& DataEntityType::UniqueId() { return Attribute("UniqueId"); }

DataEntityType& DataEntityType::UniqueId(const std::string& uuid)
{
    Attribute("UniqueId", uuid);
    return *this;
}

const std::string& DataEntityType::ValueDataType() const { return Attribute("ValueDataType"); }

std::string& DataEntityType::ValueDataType() { return Attribute("ValueDataType"); }

DataEntityType& DataEntityType::ValueDataType(const std::string& valueDataType)
{
    Attribute("ValueDataType", valueDataType);
    return *this;
}

// -----------------
// IndexedDataType
// -----------------

IndexedDataType::IndexedDataType(const std::string& metatype, const std::string& filename,
                                 const std::string& label, const XsdType& xsd)
    : InputOutputDataType(metatype, filename, label, xsd)
{
}

IndexedDataType::IndexedDataType(const std::string& metatype, const std::string& filename,
                                 const std::string& label, const FromInputXml& fromInputXml,
                                 const XsdType& xsd)
    : InputOutputDataType(metatype, filename, label, fromInputXml, xsd)
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

InputOutputDataType::InputOutputDataType(const std::string& metatype, const std::string& filename,
                                         const std::string& label, const FromInputXml& fromInputXml,
                                         const XsdType& xsd)
    : StrictEntityType(metatype, label, fromInputXml, xsd)
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
        transformedMetatype + "-" + TimeUtils::ToDataSetFormat(TimeUtils::CurrentTime());
    TimeStampedName(tsn);

    // UniqueId
    UniqueId(GenerateUuid());
}

StrictEntityType::StrictEntityType(const std::string& /*metatype*/, const std::string& label,
                                   const FromInputXml& fromInputXml, const XsdType& xsd)
    : BaseEntityType(label, fromInputXml, xsd)
{
}

const std::string& StrictEntityType::MetaType() const { return Attribute("MetaType"); }

std::string& StrictEntityType::MetaType() { return Attribute("MetaType"); }

StrictEntityType& StrictEntityType::MetaType(const std::string& metatype)
{
    Attribute("MetaType", metatype);
    return *this;
}

const std::string& StrictEntityType::TimeStampedName() const
{
    return Attribute("TimeStampedName");
}

std::string& StrictEntityType::TimeStampedName() { return Attribute("TimeStampedName"); }

StrictEntityType& StrictEntityType::TimeStampedName(const std::string& timeStampedName)
{
    Attribute("TimeStampedName", timeStampedName);
    return *this;
}

const std::string& StrictEntityType::UniqueId() const { return Attribute("UniqueId"); }

std::string& StrictEntityType::UniqueId() { return Attribute("UniqueId"); }

StrictEntityType& StrictEntityType::UniqueId(const std::string& uuid)
{
    Attribute("UniqueId", uuid);
    return *this;
}

}  // namespace internal
}  // namespace BAM
}  // namespace PacBio
