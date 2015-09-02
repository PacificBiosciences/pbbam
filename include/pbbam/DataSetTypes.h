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

#ifndef DATASETTYPES_H
#define DATASETTYPES_H

#include "pbbam/BamFile.h"
#include "pbbam/Config.h"
#include "pbbam/DataSetXsd.h"
#include "pbbam/internal/DataSetBaseTypes.h"
#include <string>

namespace PacBio {
namespace BAM {

class PBBAM_EXPORT DataSetMetadata : public internal::DataSetElement
{
public:
    DataSetMetadata(const std::string& numRecords,
                    const std::string& totalLength);

public:
    DataSetMetadata& operator+=(const DataSetMetadata& other);

public:
    const std::string& NumRecords(void) const;
    const std::string& TotalLength(void) const;
    const PacBio::BAM::Provenance& Provenance(void) const;

    std::string& NumRecords(void);
    std::string& TotalLength(void);
    PacBio::BAM::Provenance& Provenance(void);

    DataSetMetadata& NumRecords(const std::string& numRecords);
    DataSetMetadata& TotalLength(const std::string& totalLength);
    DataSetMetadata& Provenance(const PacBio::BAM::Provenance& provenance);
};

class PBBAM_EXPORT ExtensionElement : public internal::DataSetElement  {
public:
    ExtensionElement(void);
};

class PBBAM_EXPORT Extensions : public internal::DataSetListElement<ExtensionElement>
{
public:
    Extensions(void);
};

class ExternalResources;

class PBBAM_EXPORT ExternalResource : public internal::IndexedDataType
{
public:
    ExternalResource(const BamFile& bamFile);
    ExternalResource(const std::string& metatype,
                     const std::string& filename);

public:
    BamFile ToBamFile(void) const;

public:
    const PacBio::BAM::ExternalResources& ExternalResources(void) const;
    PacBio::BAM::ExternalResources& ExternalResources(void);
    ExternalResource& ExternalResources(const PacBio::BAM::ExternalResources& resources);
};

class PBBAM_EXPORT ExternalResources : public internal::DataSetListElement<ExternalResource>
{
public:
    ExternalResources(void);

    ExternalResources& operator+=(const ExternalResources& other);

public:
    void Add(const ExternalResource& ext);
    void Remove(const ExternalResource& ext);

public:
    std::vector<BamFile> BamFiles(void) const;
};

class PBBAM_EXPORT FileIndex : public internal::InputOutputDataType
{
public:
    FileIndex(const std::string& metatype, 
              const std::string& filename);
};

class PBBAM_EXPORT FileIndices : public internal::DataSetListElement<FileIndex>
{
public:
    FileIndices(void);

    void Add(const FileIndex& index);
    void Remove(const FileIndex& index);
};

class PBBAM_EXPORT Filter : public internal::DataSetElement
{
public:
    Filter(void);

public:
    const PacBio::BAM::Properties& Properties(void) const;
    PacBio::BAM::Properties& Properties(void);
    Filter& Properties(const PacBio::BAM::Properties& properties);
};

class PBBAM_EXPORT Filters : public internal::DataSetListElement<Filter>
{
public:
    Filters(void);

    Filters& operator+=(const Filters& other);

    void Add(const Filter& filter);
    void Remove(const Filter& filter);
};

class PBBAM_EXPORT ParentTool : public internal::BaseEntityType {
public:
    ParentTool(void);
};

class PBBAM_EXPORT Property : public internal::DataSetElement
{
public:
    Property(const std::string& name,
             const std::string& value,
             const std::string& op);

public:
    const std::string& Name(void) const;
    const std::string& Operator(void) const;
    const std::string& Value(void) const;

    std::string& Name(void);
    std::string& Operator(void);
    std::string& Value(void);

    Property& Name(const std::string& name);
    Property& Operator(const std::string& op);
    Property& Value(const std::string& value);
};

class PBBAM_EXPORT Properties : public internal::DataSetListElement<Property>
{
public:
    Properties(void);

    void Add(const Property& property);
    void Remove(const Property& property);
};

class PBBAM_EXPORT Provenance : public internal::DataSetElement
{
public:
    Provenance(void);

public:
    const std::string& CreatedBy(void) const;
    const std::string& CommonServicesInstanceId(void) const;
    const std::string& CreatorUserId(void) const;
    const std::string& ParentJobId(void) const;
    const PacBio::BAM::ParentTool& ParentTool(void) const;

    std::string& CreatedBy(void);
    std::string& CommonServicesInstanceId(void);
    std::string& CreatorUserId(void);
    std::string& ParentJobId(void);
    PacBio::BAM::ParentTool& ParentTool(void);

    Provenance& CreatedBy(const std::string& createdBy);
    Provenance& CommonServicesInstanceId(const std::string& id);
    Provenance& CreatorUserId(const std::string& id);
    Provenance& ParentJobId(const std::string& id);
    Provenance& ParentTool(const PacBio::BAM::ParentTool& tool);
};

class SubDataSets;

class PBBAM_EXPORT DataSetBase : public internal::StrictEntityType
{
public:
    static std::shared_ptr<DataSetBase> Create(const std::string& typeName);

public:
    DataSetBase(void);

protected:
    DataSetBase(const std::string& metatype, 
                const std::string& label, 
                const XsdType& xsd);
    DataSetBase* DeepCopy(void) const;

public:
    DataSetBase& operator+=(const DataSetBase& other);

public:
    const PacBio::BAM::ExternalResources& ExternalResources(void) const;
    const PacBio::BAM::Filters& Filters(void) const;
    const PacBio::BAM::DataSetMetadata& Metadata(void) const;
    const PacBio::BAM::SubDataSets& SubDataSets(void) const;

    PacBio::BAM::ExternalResources& ExternalResources(void);
    PacBio::BAM::Filters& Filters(void);
    PacBio::BAM::DataSetMetadata& Metadata(void);
    PacBio::BAM::SubDataSets& SubDataSets(void);

    DataSetBase& ExternalResources(const PacBio::BAM::ExternalResources& resources);
    DataSetBase& Filters(const PacBio::BAM::Filters& filters);
    DataSetBase& Metadata(const PacBio::BAM::DataSetMetadata& metadata);
    DataSetBase& SubDataSets(const PacBio::BAM::SubDataSets& subdatasets);


public:
    const NamespaceRegistry& Namespaces(void) const;
    NamespaceRegistry& Namespaces(void);

private:
    NamespaceRegistry registry_;
};

class PBBAM_EXPORT AlignmentSet : public DataSetBase
{
public:
    AlignmentSet(void);
};

class PBBAM_EXPORT BarcodeSet : public DataSetBase
{
public:
    BarcodeSet(void);
};

class PBBAM_EXPORT ConsensusAlignmentSet : public DataSetBase
{
public:
    ConsensusAlignmentSet(void);
};

class PBBAM_EXPORT ConsensusReadSet : public DataSetBase
{
public:
    ConsensusReadSet(void);
};

class PBBAM_EXPORT ContigSet : public DataSetBase
{
public:
    ContigSet(void);
};

class PBBAM_EXPORT HdfSubreadSet : public DataSetBase
{
public:
    HdfSubreadSet(void);
};

class PBBAM_EXPORT ReferenceSet : public DataSetBase
{
public:
    ReferenceSet(void);
};

class PBBAM_EXPORT SubDataSets : public internal::DataSetListElement<DataSetBase>
{
public:
    SubDataSets(void);

    SubDataSets& operator+=(const DataSetBase& other); // single
    SubDataSets& operator+=(const SubDataSets& other); // list

    void Add(const DataSetBase& subdataset);
    void Remove(const DataSetBase& subdataset);
};

class PBBAM_EXPORT SubreadSet : public DataSetBase
{
public:
    SubreadSet(void);
};

} // namespace BAM
} // namespace PacBio

#include "internal/DataSetTypes.inl"

#endif // DATASETTYPES_H
