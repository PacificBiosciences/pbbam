#ifndef PBBAM_DATASETTYPES_H
#define PBBAM_DATASETTYPES_H

#include <pbbam/Config.h>

#include <pbbam/BamFile.h>
#include <pbbam/CollectionMetadata.h>
#include <pbbam/DataSetXsd.h>
#include <pbbam/internal/DataSetBaseTypes.h>

#include <iosfwd>
#include <string>

namespace PacBio {
namespace BAM {

///
/// Set filepath resolving mode for XML output. Default is to always 'absolutize'
/// paths. Selecting 'ALLOW_RELATIVE' leaves file names verbatim.
///
enum class DataSetPathMode
{
    ABSOLUTE,
    ALLOW_RELATIVE
};

///
/// \brief The DNABarcode class represents a %DNABarcode element in
///        DataSetXML, consisting of a Name and optional UniqueId.
///
class PBBAM_EXPORT DNABarcode : public internal::DataSetElement
{
public:
    DNABarcode(const std::string& name);
    DNABarcode(const std::string& name, const std::string& uuid);
    DNABarcode(const std::string& name, const internal::FromInputXml& fromInputXml);
    DNABarcode(const std::string& name, const std::string& uuid,
               const internal::FromInputXml& fromInputXml);

    const std::string& Name() const;
    std::string& Name();
    DNABarcode& Name(const std::string& name);

    const std::string& UniqueId() const;
    std::string& UniqueId();
    DNABarcode& UniqueId(const std::string& name);
};

/// \brief The DNABarcodes class represents an %DNABarcodes element in DataSetXML.
///
/// The DNABarcodes element is essentially just a list of DNABarcode
/// objects.
///
class PBBAM_EXPORT DNABarcodes : public internal::DataSetElement
{
public:
    DNABarcodes();
    DNABarcodes(const internal::FromInputXml& fromInputXml);

public:
    void Add(const DNABarcode& barcode);
    void Remove(const DNABarcode& barcode);

public:
    using value_type = DNABarcode;
    using iterator_type = internal::DataSetElementIterator<value_type>;
    using const_iterator_type = internal::DataSetElementConstIterator<value_type>;

    const value_type& operator[](size_t index) const;
    value_type& operator[](size_t index);

    iterator_type begin();
    const_iterator_type begin() const;
    const_iterator_type cbegin() const;
    iterator_type end();
    const_iterator_type end() const;
    const_iterator_type cend() const;
};

class PBBAM_EXPORT BioSample : public internal::DataSetElement
{
public:
    BioSample(const std::string& name);
    BioSample(const std::string& name, const internal::FromInputXml& fromInputXml);

    const BAM::DNABarcodes& DNABarcodes() const;
    BAM::DNABarcodes& DNABarcodes();
    BioSample& DNABarcodes(const BAM::DNABarcodes& barcodes);

    const std::string& Name() const;
    std::string& Name();
    BioSample& Name(const std::string& name);
};

/// \brief The DNABarcodes class represents an %DNABarcodes element in DataSetXML.
///
/// The DNABarcodes element is essentially just a list of DNABarcode
/// objects.
///
class PBBAM_EXPORT BioSamples : public internal::DataSetElement
{
public:
    BioSamples();
    BioSamples(const internal::FromInputXml& fromInputXml);

public:
    void Add(const BioSample& sample);
    void Remove(const BioSample& sample);

public:
    using value_type = BioSample;
    using iterator_type = internal::DataSetElementIterator<value_type>;
    using const_iterator_type = internal::DataSetElementConstIterator<value_type>;

    const value_type& operator[](size_t index) const;
    value_type& operator[](size_t index);

    iterator_type begin();
    const_iterator_type begin() const;
    const_iterator_type cbegin() const;
    iterator_type end();
    const_iterator_type end() const;
    const_iterator_type cend() const;
};

/// \brief The ExtensionElement class represents an %ExtensionElement element in
///        DataSetXML.
///
class PBBAM_EXPORT ExtensionElement : public internal::DataSetElement
{
public:
    ExtensionElement();
    ExtensionElement(const internal::FromInputXml& fromInputXml);
};

/// \brief The Extensions class represents an %Extensions element in DataSetXML.
///
/// The Extensions element is essentially just a list of ExtensionElement
/// objects.
///
class PBBAM_EXPORT Extensions : public internal::DataSetElement
{
public:
    /// \brief Creates an empty extensions list.
    Extensions();
    Extensions(const internal::FromInputXml& fromInputXml);

public:
    using value_type = ExtensionElement;
    using iterator_type = internal::DataSetElementIterator<value_type>;
    using const_iterator_type = internal::DataSetElementConstIterator<value_type>;

    const value_type& operator[](size_t index) const;
    value_type& operator[](size_t index);

    iterator_type begin();
    const_iterator_type begin() const;
    const_iterator_type cbegin() const;
    iterator_type end();
    const_iterator_type end() const;
    const_iterator_type cend() const;
};

class ExternalResources;

/// \brief The ExternalResource class represents an %ExternalResource element in
///        DataSetXML.
///
/// An ExternalResource can itself have a child element, ExternalResources, that
/// lists related files (e.g. index files).
///
class PBBAM_EXPORT ExternalResource : public internal::IndexedDataType
{
public:
    /// \brief Creates an ExternalResource from a BamFile object.
    ///
    /// The metatype & resourceId are automatically set.
    ///
    ExternalResource(const BamFile& bamFile);

    /// \brief Creates an ExternalResource with provided \p metatype and
    ///        \p filename as resource ID.
    ///
    ExternalResource(const std::string& metatype, const std::string& filename);

    ExternalResource(const std::string& metatype, const std::string& filename,
                     const internal::FromInputXml& fromInputXml);

public:
    /// \brief Fetches the resource's ExternalResources child element.
    ///
    /// \returns const reference to child element
    /// \throws std::runtime_error if element does not exist
    ///
    const BAM::ExternalResources& ExternalResources() const;

public:
    /// \brief Fetches the resource's ExternalResources child element.
    ///
    /// This element will be created if it does not yet exist.
    ///
    /// \returns non-const reference to child element
    ///
    BAM::ExternalResources& ExternalResources();

    /// \brief Sets this resource's ExternalResources child element.
    ///
    /// This element will be created if it does not yet exist.
    ///
    /// \param[in] resources  new value for the element
    /// \returns reference to this resource object
    ///
    ExternalResource& ExternalResources(const BAM::ExternalResources& resources);

public:
    /// \brief Converts an ExternalResource to a BamFile object
    ///
    /// \returns corresponding BamFile object for this ExternalResource
    /// \throws std::runtime_error if fails to open %BAM file (e.g. does not
    ///         exist, not a %BAM file, etc.)
    ///
    /// \deprecated Use the results from DataSet::BamFiles instead. This method
    ///             cannot resolve relative filepaths and will be removed in the
    ///             near future.
    ///
    BamFile ToBamFile() const;
};

/// \brief The ExternalResources class represents an %ExternalResources element
///        in DataSetXML.
///
/// The ExternalResources element is essentially just a list of ExternalResource
/// elements.
///
class PBBAM_EXPORT ExternalResources : public internal::DataSetElement
{
public:
    /// \brief Creates an empty resource list.
    ExternalResources();
    ExternalResources(const internal::FromInputXml& fromInputXml);

    /// \brief Merges \p other resource list with this one.
    ExternalResources& operator+=(const ExternalResources& other);

public:
    /// \brief Adds an ExternalResource to this list.
    void Add(const ExternalResource& ext);

    /// \brief Removes an ExternalResource from this list.
    void Remove(const ExternalResource& ext);

public:
    /// \brief Converts resource list to BamFile objects.
    ///
    /// \deprecated Use DataSet::BamFiles instead. This method cannot resolve
    ///             relative filepaths and will be removed in the near future.
    ///
    std::vector<BamFile> BamFiles() const;

public:
    using value_type = ExternalResource;
    using iterator_type = internal::DataSetElementIterator<value_type>;
    using const_iterator_type = internal::DataSetElementConstIterator<value_type>;

    const value_type& operator[](size_t index) const;
    value_type& operator[](size_t index);

    iterator_type begin();
    const_iterator_type begin() const;
    const_iterator_type cbegin() const;
    iterator_type end();
    const_iterator_type end() const;
    const_iterator_type cend() const;
};

/// \brief The FileIndex class represents a %FileIndex element in DataSetXML.
///
/// A FileIndex is used as an auxiliary to an ExternalResource, providing
/// information about a data file's index file (e.g. for %BAM files, *.bai or
/// *.pbi).
///
class PBBAM_EXPORT FileIndex : public internal::InputOutputDataType
{
public:
    /// \brief Creates a FileIndex with provided \p metatype and \p filename as
    ///        resource ID.
    ///
    FileIndex(const std::string& metatype, const std::string& filename);

    FileIndex(const std::string& metatype, const std::string& filename,
              const internal::FromInputXml& fromInputXml);
};

/// \brief The FileIndices class represents a %FileIndices element in DataSetXML.
///
/// The FileIndices element is essentially just a list of FileIndex elements,
/// providing information about a data file's index files (e.g. for %BAM files
/// this will usually be *.bai and/or *.pbi).
///
class PBBAM_EXPORT FileIndices : public internal::DataSetElement
{
public:
    /// \brief Creates an empty index list.
    FileIndices();
    FileIndices(const internal::FromInputXml& fromInputXml);

public:
    /// \brief Adds a FileIndex to this list.
    void Add(const FileIndex& index);

    /// \brief Removes a FileIndex from this list.
    void Remove(const FileIndex& index);

public:
    using value_type = FileIndex;
    using iterator_type = internal::DataSetElementIterator<value_type>;
    using const_iterator_type = internal::DataSetElementConstIterator<value_type>;

    const value_type& operator[](size_t index) const;
    value_type& operator[](size_t index);

    iterator_type begin();
    const_iterator_type begin() const;
    const_iterator_type cbegin() const;
    iterator_type end();
    const_iterator_type end() const;
    const_iterator_type cend() const;
};

/// \brief The Filter class represents a %Filter element in DataSetXML.
///
/// The Filter element allows analysis pipelines to describe filters on data
/// that should be respected downstream, without needing to create filtered
/// intermediate files.
///
/// A filter consists of a list of Property elements, each of which must be
/// passed (logical AND) to pass the filter, e.g. property1 && property2 &&
/// property3.
///
class PBBAM_EXPORT Filter : public internal::DataSetElement
{
public:
    /// \brief Creates an empty filter.
    Filter();
    Filter(const internal::FromInputXml& fromInputXml);

public:
    /// \brief Fetches the filter's property list element.
    ///
    /// \returns const reference to child element
    /// \throws std::runtime_error if element does not exist
    ///
    const BAM::Properties& Properties() const;

public:
    /// \brief Fetches the filter's property list child element.
    ///
    /// This element will be created if it does not yet exist.
    ///
    /// \returns non-const reference to child element
    ///
    BAM::Properties& Properties();

    /// \brief Sets this filter's Properties child element.
    ///
    /// This element will be created if it does not yet exist.
    ///
    /// \param[in] properties new value for the element
    /// \returns reference to this filter object
    ///
    Filter& Properties(const BAM::Properties& properties);
};

/// \brief The Filters class represents a %Filters list element in DataSetXML.
///
/// The Filters element is essentially a list of Filter elements. For analysis
/// purpose, each filter is considered separately (logical OR) to consider which
/// data passes, e.g. filter1 || filter2 || filter3.
///
class PBBAM_EXPORT Filters : public internal::DataSetElement
{
public:
    /// \brief Creates an empty filter list.
    Filters();
    Filters(const internal::FromInputXml& fromInputXml);

    /// \brief Merges \p other filter list with this one.
    Filters& operator+=(const Filters& other);

public:
    /// \brief Adds a filter to this list.
    void Add(const Filter& filter);

    /// \brief Removes a filter from this list.
    void Remove(const Filter& filter);

public:
    using value_type = Filter;
    using iterator_type = internal::DataSetElementIterator<value_type>;
    using const_iterator_type = internal::DataSetElementConstIterator<value_type>;

    const value_type& operator[](size_t index) const;
    value_type& operator[](size_t index);

    iterator_type begin();
    const_iterator_type begin() const;
    const_iterator_type cbegin() const;
    iterator_type end();
    const_iterator_type end() const;
    const_iterator_type cend() const;
};

/// \brief The ParentTool class represents a %ParentTool element in DataSetXML.
///
class PBBAM_EXPORT ParentTool : public internal::BaseEntityType
{
public:
    /// \brief Creates an empty %ParentTool element.
    ParentTool();
    ParentTool(const internal::FromInputXml& fromInputXml);
};

/// \brief The Property class represents a %Property element in DataSetXML.
///
/// A Property is the primary building block of %DataSetXML filtering. The
/// %Property element describes a data record's property (or field), some value,
/// and a comparison operator.
///
/// For example, one could filter all %BAM records with a read accuracy at or
/// above 0.9. In C++ this could be constructed like:
/// \code{.cpp}
/// Property p("accuracy", "0.9", ">=");
/// \endcode
///
class PBBAM_EXPORT Property : public internal::DataSetElement
{
public:
    /// \brief Constructs a filter property.
    Property(const std::string& name, const std::string& value, const std::string& op);
    Property(const std::string& name, const std::string& value, const std::string& op,
             const internal::FromInputXml& fromInputXml);

public:
    /// \brief Fetches the value of property's Name attribute.
    ///
    /// \returns const reference to attribute value
    ///
    const std::string& Name() const;

    /// \brief Fetches the value of property's Operator attribute.
    ///
    /// \returns const reference to attribute value
    ///
    const std::string& Operator() const;

    /// \brief Fetches the value of property's Value attribute.
    ///
    /// \returns const reference to attribute value
    ///
    const std::string& Value() const;

public:
    /// \brief Fetches the value of property's Name attribute.
    ///
    /// \returns non-const reference to attribute value
    ///
    std::string& Name();

    /// \brief Fetches the value of property's Operator attribute.
    ///
    /// \returns non-const reference to attribute value
    ///
    std::string& Operator();

    /// \brief Fetches the value of property's Value attribute.
    ///
    /// \returns nonconst reference to attribute value
    ///
    std::string& Value();

public:
    /// \brief Sets this property's Name attribute.
    ///
    /// \param[in] name  new value for the attribute
    /// \returns reference to this property object
    ///
    Property& Name(const std::string& name);

    /// \brief Sets this property's Operator attribute.
    ///
    /// \param[in] op  new value for the attribute
    /// \returns reference to this property object
    ///
    Property& Operator(const std::string& op);

    /// \brief Sets this property's Value attribute.
    ///
    /// \param[in] value  new value for the attribute
    /// \returns reference to this property object
    ///
    Property& Value(const std::string& value);
};

/// \brief The Properties class represents a %Properties list element in
///        DataSetXML.
///
/// The Properties element is essentially a list of Property elements.
///
class PBBAM_EXPORT Properties : public internal::DataSetElement
{
public:
    /// \brief Creates an empty property list.
    Properties();
    Properties(const internal::FromInputXml& fromInputXml);

public:
    /// \brief Adds a property to this list.
    void Add(const Property& property);

    /// \brief Removes a property from this list.
    void Remove(const Property& property);

public:
    using value_type = Property;
    using iterator_type = internal::DataSetElementIterator<value_type>;
    using const_iterator_type = internal::DataSetElementConstIterator<value_type>;

    const value_type& operator[](size_t index) const;
    value_type& operator[](size_t index);

    iterator_type begin();
    const_iterator_type begin() const;
    const_iterator_type cbegin() const;
    iterator_type end();
    const_iterator_type end() const;
    const_iterator_type cend() const;
};

/// \brief The Provenance class represents a %Provenance element in DataSetXML.
///
class PBBAM_EXPORT Provenance : public internal::DataSetElement
{
public:
    /// \brief Creates a empty provenance element.
    Provenance();
    Provenance(const internal::FromInputXml& fromInputXml);

public:
    /// \brief Fetches the value of CreatedBy attribute.
    ///
    /// \returns const reference to attribute value (empty string if not
    ///          present)
    ///
    const std::string& CreatedBy() const;

    /// \brief Fetches the value of CommonServicesInstanceId attribute.
    ///
    /// \returns const reference to attribute value (empty string if not
    ///          present)
    ///
    const std::string& CommonServicesInstanceId() const;

    /// \brief Fetches the value of CreatorUserId attribute.
    ///
    /// \returns const reference to attribute value (empty string if not
    ///          present)
    ///
    const std::string& CreatorUserId() const;

    /// \brief Fetches the value of ParentJobId attribute.
    ///
    /// \returns const reference to attribute value (empty string if not
    ///          present)
    ///
    const std::string& ParentJobId() const;

    /// \brief Fetches the ParentTool child element.
    ///
    /// \returns const reference to child element
    /// \throws std::runtime_error if element does not exist
    ///
    const BAM::ParentTool& ParentTool() const;

public:
    /// \brief Fetches the value of CreatedBy attribute.
    ///
    /// This attribute will be created if it does not yet exist.
    ///
    /// \returns non-const reference to attribute value (empty string if this is
    ///          a new attribute)
    ///
    std::string& CreatedBy();

    /// \brief Fetches the value of CommonServicesInstanceId attribute.
    ///
    /// This attribute will be created if it does not yet exist.
    ///
    /// \returns non-const reference to attribute value (empty string if this is
    ///          a new attribute)
    ///
    std::string& CommonServicesInstanceId();

    /// \brief Fetches the value of CreatorUserId attribute.
    ///
    /// This attribute will be created if it does not yet exist.
    ///
    /// \returns non-const reference to attribute value (empty string if this is
    ///          a new attribute)
    ///
    std::string& CreatorUserId();

    /// \brief Fetches the value of ParentJobId attribute.
    ///
    /// This attribute will be created if it does not yet exist.
    ///
    /// \returns non-const reference to attribute value (empty string if this is
    ///          a new attribute)
    ///
    std::string& ParentJobId();

    /// \brief Fetches the ParentTool element element.
    ///
    /// This element will be created if it does not yet exist.
    ///
    /// \returns non-const reference to child element
    ///
    BAM::ParentTool& ParentTool();

public:
    /// \brief Sets the CreatedBy attribute.
    ///
    /// This attribute will be created if it does not yet exist.
    ///
    /// \param[in] createdBy  new value for the attribute
    /// \returns reference to this object
    ///
    Provenance& CreatedBy(const std::string& createdBy);

    /// \brief Sets the CommonServicesInstanceId attribute.
    ///
    /// This attribute will be created if it does not yet exist.
    ///
    /// \param[in] id  new value for the attribute
    /// \returns reference to this object
    ///
    Provenance& CommonServicesInstanceId(const std::string& id);

    /// \brief Sets the CreatorUserId attribute.
    ///
    /// This attribute will be created if it does not yet exist.
    ///
    /// \param[in] id  new value for the attribute
    /// \returns reference to this object
    ///
    Provenance& CreatorUserId(const std::string& id);

    /// \brief Sets the ParentJobId attribute.
    ///
    /// This attribute will be created if it does not yet exist.
    ///
    /// \param[in] id  new value for the attribute
    /// \returns reference to this object
    ///
    Provenance& ParentJobId(const std::string& id);

    /// \brief Sets the ParentTool child element.
    ///
    /// This element will be created if it does not yet exist.
    ///
    /// \param[in] tool  new value for the element
    /// \returns reference to this dataset object
    ///
    Provenance& ParentTool(const BAM::ParentTool& tool);
};

/// \brief The SupplementalResources class represents an %SupplementalResources element
///        in DataSetXML.
///
/// The SupplementalResources element is essentially just a list of ExternalResource
/// elements.
///
class PBBAM_EXPORT SupplementalResources : public internal::DataSetElement
{
public:
    /// \brief Creates an empty resource list.
    SupplementalResources();
    SupplementalResources(const internal::FromInputXml& fromInputXml);

    /// \brief Merges \p other resource list with this one.
    SupplementalResources& operator+=(const SupplementalResources& other);

public:
    /// \brief Adds an ExternalResource to this list.
    void Add(const ExternalResource& ext);

    /// \brief Removes an ExternalResource from this list.
    void Remove(const ExternalResource& ext);

public:
    using value_type = ExternalResource;
    using iterator_type = internal::DataSetElementIterator<value_type>;
    using const_iterator_type = internal::DataSetElementConstIterator<value_type>;

    const value_type& operator[](size_t index) const;
    value_type& operator[](size_t index);

    iterator_type begin();
    const_iterator_type begin() const;
    const_iterator_type cbegin() const;
    iterator_type end();
    const_iterator_type end() const;
    const_iterator_type cend() const;
};

/// \brief The DataSetMetadata class represents the %DataSetMetadata child
///        element in DataSetXML.
///
/// A few top-level elements are built-in, but as pbbam is not primarily a
/// DataSetXML API, most of the metadata hierarchy needs to be manually managed.
///
class PBBAM_EXPORT DataSetMetadata : public internal::DataSetElement
{
public:
    /// \name Constructors & Related Methods
    /// \{

    DataSetMetadata();
    DataSetMetadata(const internal::FromInputXml& fromInputXml);

    /// \brief Constructs a DataSetMetadata with required fields.
    DataSetMetadata(const std::string& numRecords, const std::string& totalLength);
    DataSetMetadata(const std::string& numRecords, const std::string& totalLength,
                    const internal::FromInputXml& fromInputXml);

    /// \}

public:
    /// \name Operators
    /// \{

    /// \brief Merges DataSetMetadata contents.
    ///
    /// Adds contents of \p other to this metadata object
    ///
    /// \param[in] other  some other metadata to add to this one
    /// \returns reference to this object
    ///
    DataSetMetadata& operator+=(const DataSetMetadata& other);

    /// \}

public:
    /// \name Child Elements
    /// \{

    /// \brief Fetches the text of the NumRecords element.
    ///
    /// \returns const reference to element text (empty string if not present)
    ///
    const std::string& NumRecords() const;

    /// \brief Fetches the text of the TotalLength element.
    ///
    /// \returns const reference to element text (empty string if not present)
    ///
    const std::string& TotalLength() const;

    /// \brief Fetches the Provenance element.
    ///
    /// \returns const reference to child element
    /// \throws std::runtime_error if element does not exist
    ///
    const BAM::Provenance& Provenance() const;

    /// \brief Fetches the BioSamples element.
    ///
    /// \returns const reference to child element
    /// \throws std::runtime_error if element does not exist
    ///
    const BAM::BioSamples& BioSamples() const;

    ///
    /// \brief Fetches the CollectionMetadata.
    ///
    /// \note Assumes 1 CollectionMetadata child for a given DataSetMetadata instance.
    ///
    /// \returns const reference to child element
    /// \throw std::runtime_error if element does not exist
    ///
    const BAM::CollectionMetadata& CollectionMetadata() const;

    /// \}

public:
    /// \name Child Elements
    /// \{

    /// \brief Fetches the text of the NumRecords element.
    ///
    /// This element will be created if it does not yet exist.
    ///
    /// \returns non-const reference to element text
    ///
    std::string& NumRecords();

    /// \brief Fetches the text of the TotalLength element.
    ///
    /// This element will be created if it does not yet exist.
    ///
    /// \returns non-const reference to element text
    ///
    std::string& TotalLength();

    /// \brief Fetches Provenance element.
    ///
    /// This element will be created if it does not yet exist.
    ///
    /// \returns non-const reference to child element
    ///
    BAM::Provenance& Provenance();

    /// \brief Fetches BioSamples element.
    ///
    /// This element will be created if it does not yet exist.
    ///
    /// \returns non-const reference to child element
    ///
    BAM::BioSamples& BioSamples();

    ///
    /// \brief Fetches the CollectionMetadata element.
    ///
    /// This element will be created if it does not yet exist.
    ///
    /// \note Assumes 1 CollectionMetadata child for a given DataSetMetadata instance.
    ///
    /// \return const CollectionMetadata&
    ///
    BAM::CollectionMetadata& CollectionMetadata();

    /// \}

public:
    /// \name Child Elements
    /// \{

    /// \brief Sets the text of the NumRecords element.
    ///
    /// This element will be created if it does not yet exist.
    ///
    /// \returns reference to this metadata object
    ///
    DataSetMetadata& NumRecords(const std::string& numRecords);

    /// \brief Sets the text of the TotalLength element.
    ///
    /// This element will be created if it does not yet exist.
    ///
    /// \returns reference to this metadata object
    ///
    DataSetMetadata& TotalLength(const std::string& totalLength);

    /// \brief Sets the Provenance child element.
    ///
    /// This element will be created if it does not yet exist.
    ///
    /// \returns reference to this metadata object
    ///
    DataSetMetadata& Provenance(const BAM::Provenance& provenance);

    /// \brief Sets the BioSamples child element.
    ///
    /// This element will be created if it does not yet exist.
    ///
    /// \returns reference to this metadata object
    ///
    DataSetMetadata& BioSamples(const BAM::BioSamples& samples);

    /// \brief Sets the CollectionMetadata child element.
    ///
    /// This element will be created if it does not yet exist.
    ///
    /// \returns reference to this metadata object
    ///
    DataSetMetadata& CollectionMetadata(const BAM::CollectionMetadata& metadata);

    /// \}
};

class SubDataSets;

/// \brief The DataSetBase class provides the attributes & child elements shared
///        by all dataset types.
///
/// Client code should not need to use this class directly. It should be
/// considered as more of an implementation detail and may in fact be removed
/// from public API in the future. The top-level DataSet is the recommended
/// entry point.
///
class PBBAM_EXPORT DataSetBase : public internal::StrictEntityType
{
public:
    /// \brief Creates a DataSetBase object, or one of its subclasses, from an
    ///        XML element name (e.g. SubreadSet)
    ///
    static std::shared_ptr<DataSetBase> Create(const std::string& typeName);
    static std::shared_ptr<DataSetBase> Create(const std::string& typeName,
                                               const internal::FromInputXml& fromInputXml);

public:
    /// \brief Creates an empty, generic DataSetBase.
    DataSetBase();
    DataSetBase(const internal::FromInputXml& fromInputXml);

protected:
    /// \brief Creates a DataSetBase with key values initialized.
    DataSetBase(const std::string& metatype, const std::string& label, const XsdType& xsd);
    DataSetBase(const std::string& metatype, const std::string& label,
                const internal::FromInputXml& fromInputXml, const XsdType& xsd);

    /// \brief Returns a new DataSetBase containing a deep copy of contents
    DataSetBase* DeepCopy() const;

public:
    /// \brief Merges dataset contents.
    ///
    /// Adds contents of \p other to this dataset object
    ///
    /// \param[in] other  some other dataset to add to this one
    /// \returns reference to this dataset object
    ///
    DataSetBase& operator+=(const DataSetBase& other);

public:
    /// \brief Fetches the dataset's ExternalResources element.
    ///
    /// \returns const reference to child element
    /// \throws std::runtime_error if element does not exist
    ///
    const BAM::ExternalResources& ExternalResources() const;

    /// \brief Fetches the dataset's Filters element.
    ///
    /// \returns const reference to child element
    ///
    const BAM::Filters& Filters() const;

    /// \brief Fetches the dataset's DataSetMetadata element.
    ///
    /// \returns const reference to child element
    ///
    const BAM::DataSetMetadata& Metadata() const;

    /// \brief Fetches the dataset's DataSets element.
    ///
    /// \returns const reference to child element
    ///
    const BAM::SubDataSets& SubDataSets() const;

    /// \brief Fetches the dataset's SupplementalResources element.
    ///
    /// \returns const reference to child element
    /// \throws std::runtime_error if element does not exist
    ///
    const BAM::SupplementalResources& SupplementalResources() const;

public:
    /// \brief Access this dataset's namespace info.
    ///
    /// \returns const reference to dataset's NamespaceRegistry
    ///
    const NamespaceRegistry& Namespaces() const;

public:
    /// \brief Fetches the dataset's ExternalResources element.
    ///
    /// This element will be created if it does not yet exist.
    ///
    /// \returns non-const reference to child element
    ///
    BAM::ExternalResources& ExternalResources();

    /// \brief Fetches the dataset's Filters element.
    ///
    /// This element will be created if it does not yet exist.
    ///
    /// \returns non-const reference to child element
    ///
    BAM::Filters& Filters();

    /// \brief Fetches the dataset's DataSetMetadata element.
    ///
    /// This element will be created if it does not yet exist.
    ///
    /// \returns non-const reference to child element
    ///
    BAM::DataSetMetadata& Metadata();

    /// \brief Fetches the dataset's DataSets element.
    ///
    /// This element will be created if it does not yet exist.
    ///
    /// \returns non-const reference to child element
    ///
    BAM::SubDataSets& SubDataSets();

    /// \brief Fetches the dataset's SupplementalResources element.
    ///
    /// This element will be created if it does not yet exist.
    ///
    /// \returns non-const reference to child element
    ///
    BAM::SupplementalResources& SupplementalResources();

public:
    /// \brief Sets this dataset's ExternalResources element.
    ///
    /// This element will be created if it does not yet exist.
    ///
    /// \param[in] resources  new value for the element
    /// \returns reference to this dataset object
    ///
    DataSetBase& ExternalResources(const BAM::ExternalResources& resources);

    /// \brief Sets this dataset's Filters element.
    ///
    /// This element will be created if it does not yet exist.
    ///
    /// \param[in] filters  new value for the element
    /// \returns reference to this dataset object
    ///
    DataSetBase& Filters(const BAM::Filters& filters);

    /// \brief Sets this dataset's DataSetMetadata element.
    ///
    /// This element will be created if it does not yet exist.
    ///
    /// \param[in] metadata  new value for the element
    /// \returns reference to this dataset object
    ///
    DataSetBase& Metadata(const BAM::DataSetMetadata& metadata);

    /// \brief Sets this dataset's DataSets element.
    ///
    /// This element will be created if it does not yet exist.
    ///
    /// \param[in] subdatasets  new value for the element
    /// \returns reference to this dataset object
    ///
    DataSetBase& SubDataSets(const BAM::SubDataSets& subdatasets);

    /// \brief Sets this dataset's SupplementalResources element.
    ///
    /// This element will be created if it does not yet exist.
    ///
    /// \param[in] resources  new value for the element
    /// \returns reference to this dataset object
    ///
    DataSetBase& SupplementalResources(const BAM::SupplementalResources& resources);

public:
    /// \brief Access this dataset's namespace info.
    ///
    /// \returns non-const reference to dataset's NamespaceRegistry
    ///
    NamespaceRegistry& Namespaces();

public:
    /// \brief Saves dataset XML to file.
    ///
    /// \param[in] outputFilename destination for XML contents
    /// \param[in] pathMode       print absolute paths or allow relative
    ///
    /// \throws std::runtime_error if file could be opened or if DataSet
    ///         elements could not be converted to XML
    ///
    void Save(const std::string& outputFilename,
              DataSetPathMode pathMode = DataSetPathMode::ABSOLUTE);

    /// \brief Saves dataset XML to output stream, e.g. std::cout,
    ///        std::ostringstream.
    ///
    /// \param[out] out         destination for XML contents
    /// \param[in]  pathMode    print absolute paths or allow relative
    ///
    /// \throws std::runtime_error if DataSet elements could not be converted to
    ///         XML
    ///
    void SaveToStream(std::ostream& out, DataSetPathMode pathMode = DataSetPathMode::ABSOLUTE);

public:
    ///
    /// \returns true if dataset was read from XML input
    ///
    bool FromInputXml() const;

    ///
    /// \brief Indicate that dataset was read from XML input
    ///
    void FromInputXml(bool ok);

    ///
    /// \returns (absolute) path for dataset
    ///
    const std::string& Path() const;

    ///
    /// \brief Set dataset path
    ///
    void Path(const std::string& path);

private:
    NamespaceRegistry registry_;
    std::string path_;
    bool fromInputXml_ = false;
};

/// \brief The AlignmentSet class represents an %AlignmentSet root element in
///        DataSetXML.
///
class PBBAM_EXPORT AlignmentSet : public DataSetBase
{
public:
    /// \brief Creates an empty AlignmentSet dataset.
    AlignmentSet();
    AlignmentSet(const internal::FromInputXml& fromInputXml);
};

/// \brief The BarcodeSet class represents a %BarcodeSet root element in
///        DataSetXML.
///
class PBBAM_EXPORT BarcodeSet : public DataSetBase
{
public:
    /// \brief Creates an empty BarcodeSet dataset.
    BarcodeSet();
    BarcodeSet(const internal::FromInputXml& fromInputXml);
};

/// \brief The ConsensusAlignmentSet class represents a %ConsensusAlignmentSet
///        root element in DataSetXML.
///
class PBBAM_EXPORT ConsensusAlignmentSet : public DataSetBase
{
public:
    /// \brief Creates an empty ConsensusAlignmentSet dataset.
    ConsensusAlignmentSet();
    ConsensusAlignmentSet(const internal::FromInputXml& fromInputXml);
};

/// \brief The ConsensusReadSet class represents a %ConsensusReadSet root
///        element in DataSetXML.
///
class PBBAM_EXPORT ConsensusReadSet : public DataSetBase
{
public:
    /// \brief Creates an empty ConsensusReadSet dataset.
    ConsensusReadSet();
    ConsensusReadSet(const internal::FromInputXml& fromInputXml);
};

/// \brief The ContigSet class represents a %ContigSet root element in
///        DataSetXML.
///
class PBBAM_EXPORT ContigSet : public DataSetBase
{
public:
    /// \brief Creates an empty ContigSet dataset.
    ContigSet();
    ContigSet(const internal::FromInputXml& fromInputXml);
};

/// \brief The HdfSubreadSet class represents a %HdfSubreadSet root element in
///        DataSetXML.
///
class PBBAM_EXPORT HdfSubreadSet : public DataSetBase
{
public:
    /// \brief Creates an empty HdfSubreadSet dataset.
    HdfSubreadSet();
    HdfSubreadSet(const internal::FromInputXml& fromInputXml);
};

/// \brief The ReferenceSet class represents a %ReferenceSet root element in
///        DataSetXML.
///
class PBBAM_EXPORT ReferenceSet : public DataSetBase
{
public:
    /// \brief Creates an empty ReferenceSet dataset.
    ReferenceSet();
    ReferenceSet(const internal::FromInputXml& fromInputXml);
};

/// \brief The SubDataSets class represents a %DataSets list element in
///        DataSetXML.
///
/// The SubDataSets element is essentially a list of DataSets.
///
class PBBAM_EXPORT SubDataSets : public internal::DataSetElement
{
public:
    /// \brief Creates an empty list of sub-datasets.
    SubDataSets();
    SubDataSets(const internal::FromInputXml& fromInputXml);

public:
    /// \brief Adds \p other sub-dataset to this list.
    SubDataSets& operator+=(const DataSetBase& other);  // single

    /// \brief Adds \p other sub-dataset list to this list.
    SubDataSets& operator+=(const SubDataSets& other);  // list

public:
    /// \brief Adds a sub-dataset to this list.
    void Add(const DataSetBase& subdataset);

    /// \brief Removes a sub-dataset from this list.
    void Remove(const DataSetBase& subdataset);

public:
    using value_type = DataSetBase;
    using iterator_type = internal::DataSetElementIterator<value_type>;
    using const_iterator_type = internal::DataSetElementConstIterator<value_type>;

    const value_type& operator[](size_t index) const;
    value_type& operator[](size_t index);

    iterator_type begin();
    const_iterator_type begin() const;
    const_iterator_type cbegin() const;
    iterator_type end();
    const_iterator_type end() const;
    const_iterator_type cend() const;
};

/// \brief The SubreadSet class represents a %SubreadSet root element in
///        DataSetXML.
///
class PBBAM_EXPORT SubreadSet : public DataSetBase
{
public:
    /// \brief Creates an empty SubreadSet dataset.
    SubreadSet();
    SubreadSet(const internal::FromInputXml& fromInputXml);
};

/// \brief The TranscriptSet class represents a %TranscriptSet root element in
///        DataSetXML.
///
class PBBAM_EXPORT TranscriptSet : public DataSetBase
{
public:
    /// \brief Creates an empty TranscriptSet dataset.
    TranscriptSet();
    TranscriptSet(const internal::FromInputXml& fromInputXml);
};

/// \brief The TranscriptAlignmentSet class represents a %TranscriptAlignmentSet
///        root element in DataSetXML.
///
class PBBAM_EXPORT TranscriptAlignmentSet : public DataSetBase
{
public:
    /// \brief Creates an empty TranscriptAlignmentSet dataset.
    TranscriptAlignmentSet();
    TranscriptAlignmentSet(const internal::FromInputXml& fromInputXml);
};

enum class XmlElementType
{
    GENERIC_ELEMENT,
    DATASET_METADATA,
    AUTOMATION,
    AUTOMATION_PARAMETER,
    AUTOMATION_PARAMETERS,
    BINDING_KIT,
    BIOSAMPLE,
    BIOSAMPLES,
    DNA_BARCODE,
    DNA_BARCODES,
    COLLECTIONS,
    COLLECTION_METADATA,
    CONTROL_KIT,
    EXTENSION,
    EXTENSIONS,
    EXTERNAL_RESOURCE,
    EXTERNAL_RESOURCES,
    FILE_INDEX,
    FILE_INDICES,
    FILTER,
    FILTERS,
    PARENT_TOOL,
    PPACONFIG,
    PROPERTY,
    PROPERTIES,
    PROVENANCE,
    SEQUENCING_KIT_PLATE,
    SUPPLEMENTAL_RESOURCES,
    TEMPLATE_PREP_KIT,

    GENERIC_DATASET,
    ALIGNMENT_SET,
    BARCODE_SET,
    CONSENSUS_ALIGNMENT_SET,
    CONSENSUS_READ_SET,
    CONTIG_SET,
    HDF_SUBREAD_SET,
    REFERENCE_SET,
    SUBREAD_SET,
    TRANSCRIPT_SET,
    TRANSCRIPT_ALIGNMENT_SET,
    SUBDATASETS
};

/// \returns the enum value for the requested XML element
///          (generic if not a built-in element type)
XmlElementType ElementTypeFromName(const std::string& name);

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_DATASETTYPES_H
