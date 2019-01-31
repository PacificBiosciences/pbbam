// File Description
/// \file DataSetTypes.inl
/// \brief Inline implementations for the public DataSet component classes.
//
// Author: Derek Barnett

#include "pbbam/DataSetTypes.h"

namespace PacBio {
namespace BAM {

// -------------
// DataSetBase
// -------------

inline const NamespaceRegistry& DataSetBase::Namespaces() const { return registry_; }

inline NamespaceRegistry& DataSetBase::Namespaces() { return registry_; }

// -----------------
// DataSetMetadata
// -----------------

inline const std::string& DataSetMetadata::NumRecords() const { return ChildText("NumRecords"); }

inline std::string& DataSetMetadata::NumRecords() { return ChildText("NumRecords"); }

inline DataSetMetadata& DataSetMetadata::NumRecords(const std::string& numRecords)
{
    ChildText("NumRecords", numRecords);
    return *this;
}

inline const std::string& DataSetMetadata::TotalLength() const { return ChildText("TotalLength"); }

inline std::string& DataSetMetadata::TotalLength() { return ChildText("TotalLength"); }

inline DataSetMetadata& DataSetMetadata::TotalLength(const std::string& totalLength)
{
    ChildText("TotalLength", totalLength);
    return *this;
}

// ------------
// Extensions
// ------------

inline Extensions::iterator_type Extensions::begin() { return Extensions::iterator_type(this, 0); }

inline Extensions::const_iterator_type Extensions::begin() const { return cbegin(); }

inline Extensions::const_iterator_type Extensions::cbegin() const
{
    return Extensions::const_iterator_type(this, 0);
}

inline Extensions::iterator_type Extensions::end()
{
    return Extensions::iterator_type(this, NumChildren());
}

inline Extensions::const_iterator_type Extensions::end() const { return cend(); }

inline Extensions::const_iterator_type Extensions::cend() const
{
    return Extensions::const_iterator_type(this, NumChildren());
}

inline const Extensions::value_type& Extensions::operator[](size_t index) const
{
    return dynamic_cast<const Extensions::value_type&>(*(children_.at(index).get()));
}

inline Extensions::value_type& Extensions::operator[](size_t index)
{
    return dynamic_cast<Extensions::value_type&>(*(children_.at(index).get()));
}

// -------------------
// ExternalResources
// -------------------

inline ExternalResources::iterator_type ExternalResources::begin()
{
    return ExternalResources::iterator_type(this, 0);
}

inline ExternalResources::const_iterator_type ExternalResources::begin() const { return cbegin(); }

inline ExternalResources::const_iterator_type ExternalResources::cbegin() const
{
    return ExternalResources::const_iterator_type(this, 0);
}

inline ExternalResources::iterator_type ExternalResources::end()
{
    return ExternalResources::iterator_type(this, NumChildren());
}

inline ExternalResources::const_iterator_type ExternalResources::end() const { return cend(); }

inline ExternalResources::const_iterator_type ExternalResources::cend() const
{
    return ExternalResources::const_iterator_type(this, NumChildren());
}

inline const ExternalResources::value_type& ExternalResources::operator[](size_t index) const
{
    return dynamic_cast<const ExternalResources::value_type&>(*(children_.at(index).get()));
}

inline ExternalResources::value_type& ExternalResources::operator[](size_t index)
{
    return dynamic_cast<ExternalResources::value_type&>(*(children_.at(index).get()));
}

// -------------
// FileIndices
// -------------

inline FileIndices::iterator_type FileIndices::begin()
{
    return FileIndices::iterator_type(this, 0);
}

inline FileIndices::const_iterator_type FileIndices::begin() const { return cbegin(); }

inline FileIndices::const_iterator_type FileIndices::cbegin() const
{
    return FileIndices::const_iterator_type(this, 0);
}

inline FileIndices::iterator_type FileIndices::end()
{
    return FileIndices::iterator_type(this, NumChildren());
}

inline FileIndices::const_iterator_type FileIndices::end() const { return cend(); }

inline FileIndices::const_iterator_type FileIndices::cend() const
{
    return FileIndices::const_iterator_type(this, NumChildren());
}

inline const FileIndices::value_type& FileIndices::operator[](size_t index) const
{
    return dynamic_cast<const FileIndices::value_type&>(*(children_.at(index).get()));
}

inline FileIndices::value_type& FileIndices::operator[](size_t index)
{
    return dynamic_cast<FileIndices::value_type&>(*(children_.at(index).get()));
}

// ---------
// Filters
// ---------

inline Filters::iterator_type Filters::begin() { return Filters::iterator_type(this, 0); }

inline Filters::const_iterator_type Filters::begin() const { return cbegin(); }

inline Filters::const_iterator_type Filters::cbegin() const
{
    return Filters::const_iterator_type(this, 0);
}

inline Filters::iterator_type Filters::end() { return Filters::iterator_type(this, NumChildren()); }

inline Filters::const_iterator_type Filters::end() const { return cend(); }

inline Filters::const_iterator_type Filters::cend() const
{
    return Filters::const_iterator_type(this, NumChildren());
}

inline const Filters::value_type& Filters::operator[](size_t index) const
{
    return dynamic_cast<const Filters::value_type&>(*(children_.at(index).get()));
}

inline Filters::value_type& Filters::operator[](size_t index)
{
    return dynamic_cast<Filters::value_type&>(*(children_.at(index).get()));
}

// ----------
// Property
// ----------

inline const std::string& Property::Name() const { return Attribute("Name"); }

inline std::string& Property::Name() { return Attribute("Name"); }

inline Property& Property::Name(const std::string& name)
{
    Attribute("Name", name);
    return *this;
}

inline const std::string& Property::Operator() const { return Attribute("Operator"); }

inline std::string& Property::Operator() { return Attribute("Operator"); }

inline Property& Property::Operator(const std::string& op)
{
    Attribute("Operator", op);
    return *this;
}

inline const std::string& Property::Value() const { return Attribute("Value"); }

inline std::string& Property::Value() { return Attribute("Value"); }

inline Property& Property::Value(const std::string& value)
{
    Attribute("Value", value);
    return *this;
}

// ------------
// Properties
// ------------

inline Properties::iterator_type Properties::begin() { return Properties::iterator_type(this, 0); }

inline Properties::const_iterator_type Properties::begin() const { return cbegin(); }

inline Properties::const_iterator_type Properties::cbegin() const
{
    return Properties::const_iterator_type(this, 0);
}

inline Properties::iterator_type Properties::end()
{
    return Properties::iterator_type(this, NumChildren());
}

inline Properties::const_iterator_type Properties::end() const { return cend(); }

inline Properties::const_iterator_type Properties::cend() const
{
    return Properties::const_iterator_type(this, NumChildren());
}

inline const Properties::value_type& Properties::operator[](size_t index) const
{
    return dynamic_cast<const Properties::value_type&>(*(children_.at(index).get()));
}

inline Properties::value_type& Properties::operator[](size_t index)
{
    return dynamic_cast<Properties::value_type&>(*(children_.at(index).get()));
}

// ------------
// Provenance
// ------------

inline const std::string& Provenance::CreatedBy() const { return Attribute("CreatedBy"); }

inline std::string& Provenance::CreatedBy() { return Attribute("CreatedBy"); }

inline Provenance& Provenance::CreatedBy(const std::string& createdBy)
{
    Attribute("CreatedBy", createdBy);
    return *this;
}

inline const std::string& Provenance::CommonServicesInstanceId() const
{
    return ChildText("CommonServicesInstanceId");
}

inline std::string& Provenance::CommonServicesInstanceId()
{
    return ChildText("CommonServicesInstanceId");
}

inline Provenance& Provenance::CommonServicesInstanceId(const std::string& id)
{
    ChildText("CommonServicesInstanceId", id);
    return *this;
}

inline const std::string& Provenance::CreatorUserId() const { return ChildText("CreatorUserId"); }

inline std::string& Provenance::CreatorUserId() { return ChildText("CreatorUserId"); }

inline Provenance& Provenance::CreatorUserId(const std::string& id)
{
    ChildText("CreatorUserId", id);
    return *this;
}

inline const std::string& Provenance::ParentJobId() const { return ChildText("ParentJobId"); }

inline std::string& Provenance::ParentJobId() { return ChildText("ParentJobId"); }

inline Provenance& Provenance::ParentJobId(const std::string& id)
{
    ChildText("ParentJobId", id);
    return *this;
}

inline Provenance& Provenance::ParentTool(const PacBio::BAM::ParentTool& tool)
{
    ParentTool() = tool;
    return *this;
}

// -------------
// SubDataSets
// -------------

inline SubDataSets::iterator_type SubDataSets::begin()
{
    return SubDataSets::iterator_type(this, 0);
}

inline SubDataSets::const_iterator_type SubDataSets::begin() const { return cbegin(); }

inline SubDataSets::const_iterator_type SubDataSets::cbegin() const
{
    return SubDataSets::const_iterator_type(this, 0);
}

inline SubDataSets::iterator_type SubDataSets::end()
{
    return SubDataSets::iterator_type(this, NumChildren());
}

inline SubDataSets::const_iterator_type SubDataSets::end() const { return cend(); }

inline SubDataSets::const_iterator_type SubDataSets::cend() const
{
    return SubDataSets::const_iterator_type(this, NumChildren());
}

inline const SubDataSets::value_type& SubDataSets::operator[](size_t index) const
{
    return dynamic_cast<const SubDataSets::value_type&>(*(children_.at(index).get()));
}

inline SubDataSets::value_type& SubDataSets::operator[](size_t index)
{
    return dynamic_cast<SubDataSets::value_type&>(*(children_.at(index).get()));
}

}  // namespace BAM
}  // namespace PacBio
