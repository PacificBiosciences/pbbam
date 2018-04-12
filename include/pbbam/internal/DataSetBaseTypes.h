// Author: Derek Barnett

#ifndef DATASETBASETYPES_H
#define DATASETBASETYPES_H

#include <string>
#include "pbbam/Config.h"
#include "pbbam/internal/DataSetElement.h"
#include "pbbam/internal/DataSetListElement.h"

namespace PacBio {
namespace BAM {

class DataSetMetadata;
class Extensions;
class ExternalResources;
class FileIndices;
class Filters;
class Properties;
class Provenance;

namespace internal {

class BaseEntityType : public DataSetElement
{
protected:
    BaseEntityType(const std::string& label, const XsdType& xsd = XsdType::BASE_DATA_MODEL);

public:
    const std::string& CreatedAt() const;
    const std::string& Description() const;
    const PacBio::BAM::Extensions& Extensions() const;
    const std::string& Format() const;
    const std::string& ModifiedAt() const;
    const std::string& Name() const;
    const std::string& ResourceId() const;
    const std::string& Tags() const;
    const std::string& Version() const;

    std::string& CreatedAt();
    std::string& Description();
    PacBio::BAM::Extensions& Extensions();
    std::string& Format();
    std::string& ModifiedAt();
    std::string& Name();
    std::string& ResourceId();
    std::string& Tags();
    std::string& Version();

    BaseEntityType& CreatedAt(const std::string& createdAt);
    BaseEntityType& Description(const std::string& description);
    BaseEntityType& Extensions(const PacBio::BAM::Extensions& extensions);
    BaseEntityType& Format(const std::string& format);
    BaseEntityType& ModifiedAt(const std::string& modifiedAt);
    BaseEntityType& Name(const std::string& name);
    BaseEntityType& ResourceId(const std::string& resourceId);
    BaseEntityType& Tags(const std::string& tags);
    BaseEntityType& Version(const std::string& version);
};

class DataEntityType : public BaseEntityType
{
protected:
    DataEntityType(const std::string& label, const XsdType& xsd = XsdType::BASE_DATA_MODEL);

public:
    const std::string& Checksum() const;
    const std::string& EncodedValue() const;
    const std::string& MetaType() const;
    const std::string& SimpleValue() const;
    const std::string& TimeStampedName() const;
    const std::string& UniqueId() const;
    const std::string& ValueDataType() const;

    std::string& Checksum();
    std::string& EncodedValue();
    std::string& MetaType();
    std::string& SimpleValue();
    std::string& TimeStampedName();
    std::string& UniqueId();
    std::string& ValueDataType();

    DataEntityType& Checksum(const std::string& checksum);
    DataEntityType& EncodedValue(const std::string& encodedValue);
    DataEntityType& MetaType(const std::string& metatype);
    DataEntityType& SimpleValue(const std::string& simpleValue);
    DataEntityType& TimeStampedName(const std::string& timeStampedName);
    DataEntityType& UniqueId(const std::string& uuid);
    DataEntityType& ValueDataType(const std::string& valueDataType);
};

class StrictEntityType : public BaseEntityType
{
protected:
    StrictEntityType(const std::string& metatype, const std::string& label,
                     const XsdType& xsd = XsdType::BASE_DATA_MODEL);

public:
    const std::string& MetaType() const;
    const std::string& TimeStampedName() const;
    const std::string& UniqueId() const;

    std::string& MetaType();
    std::string& TimeStampedName();
    std::string& UniqueId();

    StrictEntityType& MetaType(const std::string& metatype);
    StrictEntityType& TimeStampedName(const std::string& timeStampedName);
    StrictEntityType& UniqueId(const std::string& uuid);
};

class InputOutputDataType : public StrictEntityType
{
protected:
    InputOutputDataType(const std::string& metatype, const std::string& filename,
                        const std::string& label, const XsdType& xsd = XsdType::BASE_DATA_MODEL);
};

class IndexedDataType : public InputOutputDataType
{
protected:
    IndexedDataType(const std::string& metatype, const std::string& filename,
                    const std::string& label, const XsdType& xsd = XsdType::BASE_DATA_MODEL);

public:
    const PacBio::BAM::FileIndices& FileIndices() const;
    PacBio::BAM::FileIndices& FileIndices();
    IndexedDataType& FileIndices(const PacBio::BAM::FileIndices& indices);
};

}  // namespace internal
}  // namespace BAM
}  // namespace PacBio

#include "pbbam/internal/DataSetBaseTypes.inl"

#endif  // DATASETBASETYPES_H
