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
// --------------

inline const NamespaceRegistry& DataSetBase::Namespaces() const
{ return registry_; }

inline NamespaceRegistry& DataSetBase::Namespaces()
{ return registry_; }

// ---------------------
// DataSetMetadata
// ---------------------

inline const std::string& DataSetMetadata::NumRecords() const
{ return ChildText("NumRecords"); }

inline std::string& DataSetMetadata::NumRecords()
{ return ChildText("NumRecords"); }

inline DataSetMetadata& DataSetMetadata::NumRecords(const std::string& numRecords)
{ ChildText("NumRecords", numRecords); return *this; }

inline const std::string& DataSetMetadata::TotalLength() const
{ return ChildText("TotalLength"); }

inline std::string& DataSetMetadata::TotalLength()
{ return ChildText("TotalLength"); }

inline DataSetMetadata& DataSetMetadata::TotalLength(const std::string& totalLength)
{ ChildText("TotalLength", totalLength); return *this; }

// ----------
// Property
// ----------

inline const std::string& Property::Name() const
{ return Attribute("Name"); }

inline std::string& Property::Name()
{ return Attribute("Name"); }

inline Property& Property::Name(const std::string& name)
{ Attribute("Name", name); return *this; }

inline const std::string& Property::Operator() const
{ return Attribute("Operator"); }

inline std::string& Property::Operator()
{ return Attribute("Operator"); }

inline Property& Property::Operator(const std::string& op)
{ Attribute("Operator", op); return *this; }

inline const std::string& Property::Value() const
{ return Attribute("Value"); }

inline std::string& Property::Value()
{ return Attribute("Value"); }

inline Property& Property::Value(const std::string& value)
{ Attribute("Value", value); return *this; }

// ------------
// Provenance
// ------------

inline const std::string& Provenance::CreatedBy() const
{ return Attribute("CreatedBy"); }

inline std::string& Provenance::CreatedBy()
{ return Attribute("CreatedBy"); }

inline Provenance& Provenance::CreatedBy(const std::string& createdBy)
{ Attribute("CreatedBy", createdBy); return *this; }

inline const std::string& Provenance::CommonServicesInstanceId() const
{ return ChildText("CommonServicesInstanceId"); }

inline std::string& Provenance::CommonServicesInstanceId()
{ return ChildText("CommonServicesInstanceId"); }

inline Provenance& Provenance::CommonServicesInstanceId(const std::string& id)
{ ChildText("CommonServicesInstanceId", id); return *this; }

inline const std::string& Provenance::CreatorUserId() const
{ return ChildText("CreatorUserId"); }

inline std::string& Provenance::CreatorUserId()
{ return ChildText("CreatorUserId"); }

inline Provenance& Provenance::CreatorUserId(const std::string& id)
{ ChildText("CreatorUserId", id); return *this; }

inline const std::string& Provenance::ParentJobId() const
{ return ChildText("ParentJobId"); }

inline std::string& Provenance::ParentJobId()
{ return ChildText("ParentJobId"); }

inline Provenance& Provenance::ParentJobId(const std::string& id)
{ ChildText("ParentJobId", id); return *this; }

inline Provenance& Provenance::ParentTool(const PacBio::BAM::ParentTool& tool)
{ ParentTool() = tool; return *this; }

} // namespace BAM
} // namespace PacBio
