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

#ifndef DATASET_H
#define DATASET_H

#include "pbbam/BamFile.h"
#include "pbbam/Config.h"
#include "pbbam/DataSetTypes.h"
#include <memory>
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

class PBBAM_EXPORT DataSet
{
public:

    /// \name DataSet Types
    /// \{

    enum TypeEnum {
        GENERIC = 0
      , ALIGNMENT
      , BARCODE
      , CONSENSUS_ALIGNMENT
      , CONSENSUS_READ
      , CONTIG
      , HDF_SUBREAD
      , REFERENCE
      , SUBREAD
    };

    static DataSet::TypeEnum NameToType(const std::string& typeName);

    static std::string TypeToName(const DataSet::TypeEnum& type);
    /// \}

public:

    /// \name Constructors & Related Methods
    /// \{

    DataSet(void);
    DataSet(const DataSet::TypeEnum type);
    DataSet(const BamFile& bamFile);
    DataSet(const std::string& filename);
    DataSet(const DataSet& other);
    DataSet(DataSet&& other);
    DataSet& operator=(const DataSet& other);
    DataSet& operator=(DataSet&& other);
    ~DataSet(void);

    /// Creates a DataSet from "raw" XML data.
    static DataSet FromXml(const std::string& xml);

    /// \}

public:
    /// \name Operators
    /// \{

    DataSet& operator+=(const DataSet& other);

    /// \}

public:
    void Save(const std::string& outputFilename);
    void SaveToStream(std::ostream& out);

public:

    /// \name Attributes
    /// \{
    ///

    const std::string& Attribute(const std::string& name) const;
    std::string& Attribute(const std::string& name);
    DataSet& Attribute(const std::string& name, const std::string& value);

    const std::string& CreatedAt(void) const;
    const PacBio::BAM::Extensions& Extensions(void) const;
    const std::string& Format(void) const;
    const std::string& MetaType(void) const;
    const std::string& ModifiedAt(void) const;
    const std::string& Name(void) const;
    const std::string& ResourceId(void) const;
    const std::string& Tags(void) const;
    const std::string& TimeStampedName(void) const;
    const std::string& UniqueId(void) const;
    const std::string& Version(void) const;

    PacBio::BAM::DataSet::TypeEnum Type(void) const;
    const std::string& TypeName(void) const;

    /// \}

public:
    /// \name Child Elements
    /// \{

    const PacBio::BAM::ExternalResources& ExternalResources(void) const;
    const PacBio::BAM::Filters& Filters(void) const;
    const PacBio::BAM::DataSetMetadata& Metadata(void) const;
    const PacBio::BAM::SubDataSets& SubDataSets(void) const;

    /// \}

public:
    /// \name Attributes
    /// \{

    std::string& CreatedAt(void);
    PacBio::BAM::Extensions& Extensions(void);
    std::string& Format(void);
    std::string& MetaType(void);
    std::string& ModifiedAt(void);
    std::string& Name(void);
    std::string& ResourceId(void);
    std::string& Tags(void);
    std::string& TimeStampedName(void);
    std::string& UniqueId(void);
    std::string& Version(void);
    
    DataSet& CreatedAt(const std::string& createdAt);
    DataSet& Extensions(const PacBio::BAM::Extensions& extensions);
    DataSet& Format(const std::string& format);
    DataSet& MetaType(const std::string& metatype);
    DataSet& ModifiedAt(const std::string& modifiedAt);
    DataSet& Name(const std::string& name);
    DataSet& ResourceId(const std::string& resourceId);
    DataSet& Tags(const std::string& tags);
    DataSet& TimeStampedName(const std::string& timeStampedName);
    DataSet& UniqueId(const std::string& uuid);
    DataSet& Version(const std::string& version);

    DataSet& Type(const PacBio::BAM::DataSet::TypeEnum type);

    /// \}

public:
    /// \name Child Elements
    /// \{

    PacBio::BAM::ExternalResources& ExternalResources(void);
    PacBio::BAM::Filters& Filters(void);
    PacBio::BAM::DataSetMetadata& Metadata(void);
    PacBio::BAM::SubDataSets& SubDataSets(void);

    DataSet& ExternalResources(const PacBio::BAM::ExternalResources& resources);
    DataSet& Filters(const PacBio::BAM::Filters& filters);
    DataSet& Metadata(const PacBio::BAM::DataSetMetadata& metadata);
    DataSet& SubDataSets(const PacBio::BAM::SubDataSets& subdatasets);
    
    /// \}

private:
    std::unique_ptr<DataSetBase> d_;
};

} // namespace BAM
} // namespace PacBio

#include "pbbam/internal/DataSet.inl"

#endif // DATASET_H
