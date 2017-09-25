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
//
// File Description
/// \file DataSet.h
/// \brief Defines the DataSet class.
//
// Author: Derek Barnett

#ifndef DATASET_H
#define DATASET_H

#include "pbbam/BamFile.h"
#include "pbbam/Config.h"
#include "pbbam/DataSetTypes.h"
#include <chrono>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

/// \brief The DataSet class represents a %PacBio analyis dataset (e.g. from
///        XML).
///
/// \nosubgrouping
///
/// It provides resource paths, filters, and metadata associated with a dataset
/// under analysis.
///
class PBBAM_EXPORT DataSet
{
public:
    /// \name DataSet Type
    /// \{

    /// \brief This enum defines the currently-supported DataSet types.
    ///
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

    /// \brief Converts printable dataset type to type enum.
    ///
    /// \param[in] typeName printable dataset type
    /// \returns dataset type enum
    /// \throws std::runtime_error if \p typeName is unknown
    ///
    static DataSet::TypeEnum NameToType(const std::string& typeName);

    /// \brief Converts dataset type enum to printable name.
    ///
    /// \param[in] type dataset type enum
    /// \returns printable dataset type
    /// \throws std::runtime_error if \p type is unknown
    ///
    static std::string TypeToName(const DataSet::TypeEnum& type);

    /// \}

public:

    /// \name Constructors & Related Methods
    /// \{

    /// \brief Constructs an empty, generic DataSet.
    ///
    DataSet();

    /// \brief Constructs an empty DataSet of the type specified.
    ///
    /// \param[in] type dataset type
    /// \throws std::runtime_error if \p type is unknown
    ///
    DataSet(const DataSet::TypeEnum type);

    /// \brief Constructs a DataSet from a %BAM file.
    ///
    /// This currently defaults to a SubreadSet, with an ExternalResource
    /// pointing to BamFile::Filename.
    ///
    /// \param[in] bamFile  BamFile object
    ///
    DataSet(const BamFile& bamFile);

    /// \brief Loads a DataSet from a file.
    ///
    /// \p filename may be one of the following types, indicated by its extension:\n
    ///  - %BAM ("*.bam") \n
    ///  - FOFN ("*.fofn") \n
    ///  - FASTA ("*.fa" or "*.fasta") \n
    ///  - DataSetXML ("*.xml") \n
    ///
    /// \param[in] filename  input filename
    /// \throws std::runtime_error if \p filename has an unsupported extension,
    ///         or if a valid DataSet could not be created from its contents
    ///
    DataSet(const std::string& filename);

    /// \brief Constructs a DataSet from a list of files.
    ///
    /// \param[in] filenames  input filenames
    /// \throws std::runtime_error if DataSet could not be created from
    ///         \p filenames
    ///
    DataSet(const std::vector<std::string>& filenames);

    DataSet(const DataSet& other);
    DataSet(DataSet&&) = default;
    DataSet& operator=(const DataSet& other);
    DataSet& operator=(DataSet&&) = default;
    ~DataSet() = default;

    /// \brief Creates a DataSet from "raw" XML data.
    ///
    /// \param[in] xml DataSetXML text
    ///
    static DataSet FromXml(const std::string& xml);

    /// \}

public:
    /// \name Operators
    /// \{

    /// \brief Merges DataSet contents.
    ///
    /// Adds contents of \p other to this dataset object
    ///
    /// \param[in] other  some other dataset to add to this one
    /// \returns reference to this dataset object
    ///
    DataSet& operator+=(const DataSet& other);

    /// \}

public:
    /// \name Serialization
    /// \{

    /// \brief Saves dataset XML to file.
    ///
    /// \param[in] outputFilename destination for XML contents
    ///
    /// \throws std::runtime_error if file could be opened or if DataSet
    ///         elements could not be converted to XML
    ///
    void Save(const std::string& outputFilename);

    /// \brief Saves dataset XML to output stream, e.g. std::cout,
    ///        std::stringstream.
    ///
    /// \param[out] out destination for XML contents
    ///
    /// \throws std::runtime_error if DataSet elements could not be converted to
    ///         XML
    ///
    void SaveToStream(std::ostream& out);

    /// \}

public:

    /// \name Attributes
    /// \{
    ///

    /// \brief Fetches the value of a DataSet root element's attribute.
    ///
    /// These are the attributes attached to the root dataset element: \n
    /// \verbatim <SubreadSet foo="x" bar="y" /> \endverbatim
    ///
    /// Built-in accessors exist for the standard attributes (e.g. CreatedAt)
    /// but additional attributes can be used as well via these generic
    /// Attribute methods.
    ///
    /// \param[in] name root element's attribute name
    /// \returns const reference to attribute's value (empty string if not
    ///          present)
    ///
    const std::string& Attribute(const std::string& name) const;

    /// \brief Fetches the value of dataset's CreatedAt attribute.
    ///
    /// \returns const reference to attribute's value (empty string if not
    ///          present)
    ///
    const std::string& CreatedAt() const;

    /// \brief Fetches the value of dataset's Format attribute.
    ///
    /// \returns const reference to attribute's value (empty string if not
    ///          present)
    ///
    const std::string& Format() const;

    /// \brief Fetches the value of dataset's MetaType attribute.
    ///
    /// \returns const reference to attribute's value (empty string if not
    ///          present)
    ///
    const std::string& MetaType() const;

    /// \brief Fetches the value of dataset's ModifiedAt attribute.
    ///
    /// \returns const reference to attribute's value (empty string if not
    ///          present)
    ///
    const std::string& ModifiedAt() const;

    /// \brief Fetches the value of dataset's Name attribute.
    ///
    /// \returns const reference to attribute's value (empty string if not
    ///          present)
    ///
    const std::string& Name() const;

    /// \brief Fetches the value of dataset's ResourceId attribute.
    ///
    /// \returns const reference to attribute's value (empty string if not
    ///          present)
    ///
    const std::string& ResourceId() const;

    /// \brief Fetches the value of dataset's Tags attribute.
    ///
    /// \returns const reference to attribute's value (empty string if not
    ///          present)
    ///
    const std::string& Tags() const;

    /// \brief Fetches the value of dataset's TimeStampedName attribute.
    ///
    /// \returns const reference to attribute's value (empty string if not
    ///          present)
    ///
    const std::string& TimeStampedName() const;

    /// \brief Fetches the value of dataset's UniqueId attribute.
    ///
    /// \returns const reference to attribute's value (empty string if not
    ///          present)
    ///
    const std::string& UniqueId() const;

    /// \brief Fetches the value of dataset's Version attribute.
    ///
    /// \returns const reference to attribute's value (empty string if not
    ///          present)
    ///
    const std::string& Version() const;

    /// \}

public:
    /// \name DataSet Type
    /// \{

    /// \brief Fetches the dataset's type.
    ///
    /// \returns dataset type enum
    ///
    PacBio::BAM::DataSet::TypeEnum Type() const;

    /// \brief Fetches the dataset's type.
    ///
    /// \returns printable dataset type
    ///
    std::string TypeName() const;

    /// \}

public:
    /// \name Child Elements
    /// \{

    /// \brief Fetches the dataset's Extensions element.
    ///
    /// \returns const reference to child element
    /// \throws std::runtime_error if element does not exist
    ///
    const PacBio::BAM::Extensions& Extensions() const;

    /// \brief Fetches the dataset's ExternalResources element.
    ///
    /// \returns const reference to child element
    /// \throws std::runtime_error if element does not exist
    ///
    const PacBio::BAM::ExternalResources& ExternalResources() const;

    /// \brief Fetches the dataset's Filters element.
    ///
    /// \returns const reference to child element
    ///
    const PacBio::BAM::Filters& Filters() const;

    /// \brief Fetches the dataset's DataSetMetadata element.
    ///
    /// \returns const reference to child element
    ///
    const PacBio::BAM::DataSetMetadata& Metadata() const;

    /// \brief Fetches the dataset's DataSets element.
    ///
    /// \returns const reference to child element
    ///
    const PacBio::BAM::SubDataSets& SubDataSets() const;

    /// \}

public:
    /// \name Resource Handling
    /// \{
 
    /// \brief Returns all of this dataset's resource files, with relative
    ///        filepaths already resolved.
    ///
    /// Includes both primary resources (e.g. subread BAM files), as well as all 
    /// secondary or child resources (e.g. index files, scraps BAM, etc). 
    ///
    /// \returns vector of (resolveD) filepaths
    ///
    /// \sa DataSet::ResolvedResourceIds
    ///
    std::vector<std::string> AllFiles() const;
 
    /// \brief Returns this dataset's primary %BAM resources, with relative
    ///        filepaths already resolved.
    ///
    /// Primary resources are those listed as top-level %ExternalResources, not
    /// associated files (indices, references, scraps %BAMs, etc.).
    ///
    /// \returns vector of BamFiles
    ///
    /// \sa DataSet::ResolvedResourceIds
    ///
    std::vector<BamFile> BamFiles() const;

    /// \brief Returns this dataset's primary FASTA resources, with relative
    ///        filepaths already resolved.
    ///
    /// Primary resources are those listed as top-level %ExternalResources, not
    /// associated files (indices, references, scraps %BAMs, etc.).
    ///
    /// \returns vector of filepaths to FASTA resources
    ///
    /// \sa DataSet::ResolvedResourceIds
    ///
    std::vector<std::string> FastaFiles() const;

    /// \brief Returns all primary external resource filepaths, with relative
    ///        paths resolved.
    ///
    /// Primary resources are those listed as top-level %ExternalResources, not
    /// associated files (indices, references, scraps %BAMs, etc.).
    ///
    /// \sa ResolvePath
    ///
    /// \returns resourceIds
    ///
    std::vector<std::string> ResolvedResourceIds() const;

    /// \brief Resolves a filepath (that may be relative to the dataset).
    ///
    /// A DataSet's resources may be described using absolute filepaths or with
    /// relative paths. For absolute paths, nothing is changed from the input.
    /// For relative paths, these are resolved using the DataSet's own path
    /// as a starting point. A DataSet's own path will be one of:\n
    ///  1 - the location of its XML or %BAM input file, e.g. created using
    ///      DataSet("foo.xml") or DataSet("foo.bam")\n
    ///  2 - application's current working directory for all other DataSet
    ///      construction methods { DataSet(), DataSet(type),
    ///      DataSet("foo.fofn") }\n
    ///
    /// \param[in] originalPath     input file path (absolute or relative)
    /// \returns resolved path
    ///
    std::string ResolvePath(const std::string& originalPath) const;

    /// \returns sequence chemistry info for all read groups in this dataset
    ///
    /// \sa ReadGroupInfo::SequencingChemistry
    ///
    std::set<std::string> SequencingChemistries() const;

    /// \}

public:
    /// \name XML Namespace Handling
    /// \{

    /// \brief Access this dataset's namespace info.
    ///
    /// \returns const reference to dataset's NamespaceRegistry
    ///
    const NamespaceRegistry& Namespaces() const;

    /// \}

public:
    /// \name Attributes
    /// \{

    /// \brief Fetches the value of a DataSet root element's attribute.
    ///
    /// These are the attributes attached to the root dataset element: \n
    /// \verbatim <SubreadSet foo="x" bar="y" /> \endverbatim
    ///
    /// Built-in accessors exist for the standard attributes (e.g. CreatedAt)
    /// but additional attributes can be used as well via these generic methods.
    ///
    /// A new attribute will be created if it does not yet exist.
    ///
    /// \param[in] name root element's attribute name
    /// \returns non-const reference to attribute's value (empty string if this
    ///          is a new attribute)
    ///
    std::string& Attribute(const std::string& name);

    /// \brief Fetches the value of dataset's CreatedAt attribute.
    ///
    /// This attribute will be created if it does not yet exist.
    ///
    /// \returns non-const reference to attribute's value (empty string if this
    ///          is a new attribute)
    ///
    std::string& CreatedAt();

    /// \brief Fetches the value of dataset's Format attribute.
    ///
    /// This attribute will be created if it does not yet exist.
    ///
    /// \returns non-const reference to attribute's value (empty string if this
    ///          is a new attribute)
    ///
    std::string& Format();

    /// \brief Fetches the value of dataset's MetaType attribute.
    ///
    /// This attribute will be created if it does not yet exist.
    ///
    /// \returns non-const reference to attribute's value (empty string if this
    ///          is a new attribute)
    ///
    std::string& MetaType();

    /// \brief Fetches the value of dataset's ModifiedAt attribute.
    ///
    /// This attribute will be created if it does not yet exist.
    ///
    /// \returns non-const reference to attribute's value (empty string if this
    ///          is a new attribute)
    ///
    std::string& ModifiedAt();

    /// \brief Fetches the value of dataset's Name attribute.
    ///
    /// This attribute will be created if it does not yet exist.
    ///
    /// \returns non-const reference to attribute's value (empty string if this
    ///          is a new attribute)
    ///
    std::string& Name();

    /// \brief Fetches the value of dataset's ResourceId attribute.
    ///
    /// This attribute will be created if it does not yet exist.
    ///
    /// \returns non-const reference to attribute's value (empty string if this
    ///          is a new attribute)
    ///
    std::string& ResourceId();

    /// \brief Fetches the value of dataset's Tags attribute.
    ///
    /// This attribute will be created if it does not yet exist.
    ///
    /// \returns non-const reference to attribute's value (empty string if this
    ///          is a new attribute)
    ///
    std::string& Tags();

    /// \brief Fetches the value of dataset's TimeStampedName attribute.
    ///
    /// This attribute will be created if it does not yet exist.
    ///
    /// \returns non-const reference to attribute's value (empty string if this
    ///          is a new attribute)
    ///
    std::string& TimeStampedName();

    /// \brief Fetches the value of dataset's UniqueId attribute.
    ///
    /// This attribute will be created if it does not yet exist.
    ///
    /// \returns non-const reference to attribute's value (empty string if this
    ///          is a new attribute)
    ///
    std::string& UniqueId();

    /// \brief Fetches the value of dataset's Version attribute.
    ///
    /// This attribute will be created if it does not yet exist.
    ///
    /// \returns non-const reference to attribute's value (empty string if this
    ///          is a new attribute)
    ///
    std::string& Version();
    
    /// \}

public:
    /// \name Attributes
    /// \{

    /// \brief Sets this dataset's XML attribute \p name, with \p value
    ///
    /// These are the attributes attached to the root dataset element: \n
    /// \verbatim <SubreadSet foo="x" bar="y" /> \endverbatim
    ///
    /// Built-in accessors exist for the standard attributes (e.g. CreatedAt)
    /// but additional attributes can be used as well via these generic methods.
    ///
    /// The attribute will be created if it does not yet exist.
    ///
    /// \param[in] name   root element's attribute name
    /// \param[in] value  new value for the attribute
    /// \returns reference to this dataset object
    ///
    DataSet& Attribute(const std::string& name, const std::string& value);

    /// \brief Sets this dataset's CreatedAt attribute.
    ///
    /// This attribute will be created if it does not yet exist.
    ///
    /// \param[in] createdAt  new value for the attribute
    /// \returns reference to this dataset object
    ///
    DataSet& CreatedAt(const std::string& createdAt);

    /// \brief Sets this dataset's Format attribute.
    ///
    /// This attribute will be created if it does not yet exist.
    ///
    /// \param[in] format  new value for the attribute
    /// \returns reference to this dataset object
    ///
    DataSet& Format(const std::string& format);

    /// \brief Sets this dataset's MetaType attribute.
    ///
    /// This attribute will be created if it does not yet exist.
    ///
    /// \param[in] metatype  new value for the attribute
    /// \returns reference to this dataset object
    ///
    DataSet& MetaType(const std::string& metatype);

    /// \brief Sets this dataset's ModifiedAt attribute.
    ///
    /// This attribute will be created if it does not yet exist.
    ///
    /// \param[in] modifiedAt  new value for the attribute
    /// \returns reference to this dataset object
    ///
    DataSet& ModifiedAt(const std::string& modifiedAt);

    /// \brief Sets this dataset's Name attribute.
    ///
    /// This attribute will be created if it does not yet exist.
    ///
    /// \param[in] name  new value for the attribute
    /// \returns reference to this dataset object
    ///
    DataSet& Name(const std::string& name);

    /// \brief Sets this dataset's ResourceId attribute.
    ///
    /// This attribute will be created if it does not yet exist.
    ///
    /// \param[in] resourceId  new value for the attribute
    /// \returns reference to this dataset object
    ///
    DataSet& ResourceId(const std::string& resourceId);

    /// \brief Sets this dataset's Tags attribute.
    ///
    /// This attribute will be created if it does not yet exist.
    ///
    /// \param[in] tags  new value for the attribute
    /// \returns reference to this dataset object
    ///
    DataSet& Tags(const std::string& tags);

    /// \brief Sets this dataset's TimeStampedName attribute.
    ///
    /// This attribute will be created if it does not yet exist.
    ///
    /// \param[in] timeStampedName  new value for the attribute
    /// \returns reference to this dataset object
    ///
    DataSet& TimeStampedName(const std::string& timeStampedName);

    /// \brief Sets this dataset's UniqueId attribute.
    ///
    /// This attribute will be created if it does not yet exist.
    ///
    /// \param[in] uuid  new value for the attribute
    /// \returns reference to this dataset object
    ///
    DataSet& UniqueId(const std::string& uuid);

    /// \brief Sets this dataset's Version attribute.
    ///
    /// This attribute will be created if it does not yet exist.
    ///
    /// \param[in] version  new value for the attribute
    /// \returns reference to this dataset object
    ///
    DataSet& Version(const std::string& version);

    /// \}

public:
    /// \name DataSet Type
    /// \{

    /// \brief Edits dataset type.
    ///
    /// \param[in] type  new dataset type
    /// \returns reference to this dataset object
    ///
    DataSet& Type(const PacBio::BAM::DataSet::TypeEnum type);

    /// \}

public:
    /// \name Child Elements
    /// \{

    /// \brief Fetches the dataset's Extensions element.
    ///
    /// This element will be created if it does not yet exist.
    ///
    /// \returns non-const reference to child element
    ///
    PacBio::BAM::Extensions& Extensions();

    /// \brief Fetches the dataset's ExternalResources element.
    ///
    /// This element will be created if it does not yet exist.
    ///
    /// \returns non-const reference to child element
    ///
    PacBio::BAM::ExternalResources& ExternalResources();

    /// \brief Fetches the dataset's Filters element.
    ///
    /// This element will be created if it does not yet exist.
    ///
    /// \returns non-const reference to child element
    ///
    PacBio::BAM::Filters& Filters();

    /// \brief Fetches the dataset's DataSetMetadata element.
    ///
    /// This element will be created if it does not yet exist.
    ///
    /// \returns non-const reference to child element
    ///
    PacBio::BAM::DataSetMetadata& Metadata();

    /// \brief Fetches the dataset's DataSets element.
    ///
    /// This element will be created if it does not yet exist.
    ///
    /// \returns non-const reference to child element
    ///
    PacBio::BAM::SubDataSets& SubDataSets();

    /// \}

public:
    /// \name Child Elements
    /// \{

    /// \brief Sets this dataset's Extensions element.
    ///
    /// This element will be created if it does not yet exist.
    ///
    /// \param[in] extensions  new value for the element
    /// \returns reference to this dataset object
    ///
    DataSet& Extensions(const PacBio::BAM::Extensions& extensions);

    /// \brief Sets this dataset's ExternalResources element.
    ///
    /// This element will be created if it does not yet exist.
    ///
    /// \param[in] resources  new value for the element
    /// \returns reference to this dataset object
    ///
    DataSet& ExternalResources(const PacBio::BAM::ExternalResources& resources);

    /// \brief Sets this dataset's Filters element.
    ///
    /// This element will be created if it does not yet exist.
    ///
    /// \param[in] filters  new value for the element
    /// \returns reference to this dataset object
    ///
    DataSet& Filters(const PacBio::BAM::Filters& filters);

    /// \brief Sets this dataset's DataSetMetadata element.
    ///
    /// This element will be created if it does not yet exist.
    ///
    /// \param[in] metadata  new value for the element
    /// \returns reference to this dataset object
    ///
    DataSet& Metadata(const PacBio::BAM::DataSetMetadata& metadata);

    /// \brief Sets this dataset's DataSets element.
    ///
    /// This element will be created if it does not yet exist.
    ///
    /// \param[in] subdatasets  new value for the element
    /// \returns reference to this dataset object
    ///
    DataSet& SubDataSets(const PacBio::BAM::SubDataSets& subdatasets);
    
    /// \}

public:
    /// \name XML Namespace Handling
    /// \{

    /// \brief Access this dataset's namespace info.
    ///
    /// \returns non-const reference to dataset's NamespaceRegistry
    ///
    NamespaceRegistry& Namespaces();

    /// \}

private:
    std::unique_ptr<DataSetBase> d_;
    std::string path_;
};

/// \name DataSet Timestamp Utilities
/// \{

/// \brief Fetches current time, in "DataSetXML format".
///
/// \returns DataSetXML formatted timestamp
///
/// \sa ToDataSetFormat
///
PBBAM_EXPORT std::string CurrentTimestamp();

/// \brief Converts a time_point to "DataSetXML-formatted" timestamp.
///
/// This is the format used as a component of the DataSet::TimeStampedName
/// (yymmdd_HHmmssttt>.
///
/// \returns "DataSetXML-formatted" timestamp
///
PBBAM_EXPORT std::string ToDataSetFormat(const std::chrono::system_clock::time_point& tp);

/// \brief Converts a time_t to "DataSetXML-formatted" timestamp.
///
/// This is the format used as a component of the DataSet::TimeStampedName
/// (yymmdd_HHmmssttt>.
///
/// \returns "DataSetXML-formatted" timestamp
///
PBBAM_EXPORT std::string ToDataSetFormat(const time_t& tp);

/// \brief Converts a time_point to ISO-8601 formatted timestamp.
///
/// This is the format used in DataSet::CreatedAt and DataSet::ModifiedAt.
///
/// \returns ISO-8601 formatted timestamp
///
PBBAM_EXPORT std::string ToIso8601(const std::chrono::system_clock::time_point& tp);

/// \brief Converts a time_t to ISO-8601 formatted timestamp.
///
/// This is the format used in DataSet::CreatedAt and DataSet::ModifiedAt.
///
/// \returns ISO-8601 formatted timestamp
///
PBBAM_EXPORT std::string ToIso8601(const time_t& t);

/// \}

} // namespace BAM
} // namespace PacBio

#include "pbbam/internal/DataSet.inl"

#endif // DATASET_H
