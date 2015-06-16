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

#ifndef DATASETBASE_H
#define DATASETBASE_H

#include "pbbam/Config.h"
#include "pbbam/dataset/DataSetMetadata.h"
#include "pbbam/dataset/ExternalDataReferences.h"
#include "pbbam/dataset/Filters.h"
#include "pbbam/dataset/SubDataSets.h"
#include "pbbam/internal/DataSetElement.h"
#include <ostream>
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

enum class DataSetType
{
    GENERIC       = 0
  , ALIGNMENTSET
  , BARCODESET
  , CCSREADSET
  , CONTIGSET
  , REFERENCESET  = 5
  , SUBREADSET
};

class PBBAM_EXPORT DataSetBase : public internal::DataSetElement
{
public:
    /// \name DataSet Type
    /// \{

    /// \returns enum value for name string
    ///
    static DataSetType TypeForName(const std::string& name);

    /// \returns name string for enum value
    ///
    static std::string NameForType(const DataSetType type);

    /// \}

public:
    /// \name Constructors & Related Methods
    /// \{
    DataSetBase(void);

    /// Constructs an empty dataset with the specified type.
    ///
    DataSetBase(const DataSetType& type);

    /// Constructs a dataset using the contents of \p filename. This file
    /// can be either a "direct" data file (e.g. BAM) or a dataset description
    /// (e.g. DataSetXML, FOFN)
    ///
    DataSetBase(const std::string& filename);

    /// Constructs a dataset using the contents of the \p uris provided.
    /// These files can contain "direct" data (e.g. BAM) or contain
    /// dataset descriptions (e.g. DataSetXML, FOFN).
    ///
    DataSetBase(const std::vector<std::string>& uris);

    DataSetBase(const DataSetBase& other);
    DataSetBase(DataSetBase&& other);
    DataSetBase& operator=(const DataSetBase& other);
    DataSetBase& operator=(DataSetBase&& other);
    ~DataSetBase(void);

    /// \}



public:
    /// \name Header Attributes
    /// \{

    /// \returns "CreatedAt" attribute value (or empty string if none present)
    ///
    const std::string& CreatedAt(void) const;

    /// \returns "MetaType" attribute value (or empty string if none present)
    ///
    const std::string& MetaType(void) const;

    /// \returns "Name" attribute value (or empty string if none present)
    ///
    const std::string& Name(void) const;

    /// \returns "Tags" attribute value (or empty string if none present)
    ///
    const std::string& Tags(void) const;

    /// \returns enum describing dataset type
    ///
    DataSetType Type(void) const;

    /// \returns "UniqueId" attribute value (or empty string if none present)
    ///
    const std::string& UniqueId(void) const;

    /// \returns "Version" attribute value (or empty string if none present)
    ///
    const std::string& Version(void) const;

    /// \}

public:
    /// \name DataSet Components
    /// \{

    /// \returns dataset's external data reference list
    ///
    const ExternalDataReferences& ExternalDataReferenceList(void) const;

    /// \returns dataset's filter list
    ///
    const Filters& FilterList(void) const;

    /// \returns dataset's subdataset list
    ///
    const SubDataSets& SubDataSetList(void) const;

    /// \returns number of external data references
    ///
    size_t NumExternalDataReferences(void) const;

    /// \returns number of filters
    ///
    size_t NumFilters(void) const;

    /// \returns number of subdatasets
    ///
    size_t NumSubDataSets(void) const;

    /// \}

public:
    /// \n I/O Methods
    /// \{


    /// Writes DataSetXML to \p fn
    ///
    void Write(const std::string& fn) const;

    /// Writes DataSetXML to stderr
    ///
    void WriteToStderr(void) const;

    /// Writes DataSetXML to stdout
    ///
    void WriteToStdout(void) const;

    /// Writes DataSetXML to provided out stream
    ///
    void WriteToStream(std::ostream& out) const;

    /// \}

public:
    /// \name Header Attributes
    /// \{

    /// Sets this dataset's "CreatedAt" attribute.
    ///
    /// \todo: format validation
    ///
    /// \param[in] timestamp
    /// \returns reference to this dataset
    ///
    DataSetBase& CreatedAt(const std::string& timestamp);

    /// Sets this dataset's "MetaType" attribute.
    ///
    /// \param[in] metatype name
    /// \returns reference to this dataset
    ///
    DataSetBase& MetaType(const std::string& metatype);

    /// Sets this dataset's "Name" attribute.
    ///
    /// \param[in] name
    /// \returns reference to this dataset
    ///
    DataSetBase& Name(const std::string& name);

    /// Sets this dataset's "Tags" attribute.
    ///
    /// \param[in] tags string
    /// \returns reference to this dataset
    ///
    DataSetBase& Tags(const std::string& tags);

    /// Sets this dataset's type.
    ///
    /// \param[in] type
    /// \returns reference to this dataset
    ///
    DataSetBase& Type(DataSetType type);

    /// Sets this dataset's "UniqueId" attribute.
    ///
    /// \todo UUID-handling/validation
    ///
    /// \param[in] uuid
    /// \returns reference to this dataset
    ///
    DataSetBase& UniqueId(const std::string& uuid);

    /// Sets this dataset's "Version" attribute.
    ///
    /// \param[in] version
    /// \returns reference to this dataset
    ///
    DataSetBase& Version(const std::string& version);

    /// \}

public:
    /// \name DataSet Components
    /// \{

    /// Adds \p ref to external data reference list
    ///
    /// \param[in] ref
    /// \returns reference to this dataset
    ///
    DataSetBase& AddExternalDataReference(const ExternalDataReference& ref);

    /// Adds \p filter to filter list
    ///
    /// \param[in] filter
    /// \returns reference to this dataset
    ///
    DataSetBase& AddFilter(const Filter& filter);

    /// Adds \p subdataset to subdataset list
    ///
    /// \param[in] subdataset
    /// \returns reference to this dataset
    ///
    DataSetBase& AddSubDataSet(const SubDataSet& subdataset);

    /// \returns editable external data reference list
    ExternalDataReferences& ExternalDataReferenceList(void);

    /// \returns editable filters list
    Filters& FilterList(void);

    /// Removes \p ref from external data reference list
    ///
    /// \param[in] ref
    /// \returns reference to this dataset
    ///
    DataSetBase& RemoveExternalDataReference(const ExternalDataReference& ref);

    /// Removes \p filter from filter list
    ///
    /// \param[in] filter
    /// \returns reference to this dataset
    ///
    DataSetBase& RemoveFilter(const Filter& filter);

    /// Removes \p subdataset from subdataset list
    ///
    /// \param[in] subdataset
    /// \returns reference to this dataset
    ///
    DataSetBase& RemoveSubDataSet(const SubDataSet& subdataset);

    /// \returns editable sub-dataset list
    SubDataSets& SubDataSetList(void);

    /// \}

public:
    /// \name Merging Datasets
    /// \{

    /// Merges 2 datasets.
    ///
    /// Performs a union of external data references.
    ///
    /// Currently, throws exception when filters and/or metadata differ.
    /// These will likely combined more intelligently in the future.
    ///
    /// \returns merged dataset
    ///
    DataSetBase operator+(const DataSetBase& other) const;

    /// Merges \p other dataset with this dataset.
    ///
    /// Performs a union of external data references.
    ///
    /// Currently, throws exception when filters and/or metadata differ.
    /// These will likely combined more intelligently in the future.
    ///
    /// \returns reference to this dataset
    ///
    DataSetBase& operator+=(const DataSetBase& other);

    /// \}
};

} // namespace BAM
} // namespace PacBio

#endif // DATASETBASE_H
