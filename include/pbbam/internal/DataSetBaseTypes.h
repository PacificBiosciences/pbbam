// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

// Author: Derek Barnett

#ifndef DATASETBASETYPES_H
#define DATASETBASETYPES_H

#include "pbbam/Config.h"
#include "pbbam/internal/DataSetElement.h"
#include "pbbam/internal/DataSetListElement.h"
#include <string>

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
    BaseEntityType(const std::string& label,
                   const XsdType& xsd = XsdType::BASE_DATA_MODEL);

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
    DataEntityType(const std::string& label,
                   const XsdType& xsd = XsdType::BASE_DATA_MODEL);

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
    StrictEntityType(const std::string& metatype,
                     const std::string& label,
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
    InputOutputDataType(const std::string& metatype, 
                        const std::string& filename,
                        const std::string& label,
                        const XsdType& xsd = XsdType::BASE_DATA_MODEL);
};

class IndexedDataType : public InputOutputDataType
{
protected:
    IndexedDataType(const std::string& metatype, 
                    const std::string& filename,
                    const std::string& label,
                    const XsdType& xsd = XsdType::BASE_DATA_MODEL);

public:
    const PacBio::BAM::FileIndices& FileIndices() const;
    PacBio::BAM::FileIndices& FileIndices();
    IndexedDataType& FileIndices(const PacBio::BAM::FileIndices& indices);
};

} // namespace internal
} // namespace BAM
} // namespace PacBio

#include "pbbam/internal/DataSetBaseTypes.inl"

#endif // DATASETBASETYPES_H
