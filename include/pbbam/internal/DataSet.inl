// File Description
/// \file DataSet.inl
/// \brief Inline implementations for the DataSet class.
//
// Author: Derek Barnett

#include "pbbam/DataSet.h"

namespace PacBio {
namespace BAM {

inline const std::string& DataSet::Attribute(const std::string& name) const
{ return d_->Attribute(name); }

inline std::string& DataSet::Attribute(const std::string& name)
{ return d_->Attribute(name); }

inline DataSet& DataSet::Attribute(const std::string& name, const std::string& value)
{ d_->Attribute(name, value); return *this; }

inline const std::string& DataSet::CreatedAt() const
{ return d_->CreatedAt(); }

inline std::string& DataSet::CreatedAt()
{ return d_->CreatedAt(); }

inline DataSet& DataSet::CreatedAt(const std::string& createdAt)
{ d_->CreatedAt(createdAt); return *this; }

inline const PacBio::BAM::Extensions& DataSet::Extensions() const
{ return d_->Extensions(); }

inline PacBio::BAM::Extensions& DataSet::Extensions()
{ return d_->Extensions(); }

inline DataSet& DataSet::Extensions(const PacBio::BAM::Extensions& extensions)
{ d_->Extensions(extensions); return *this; }

inline const PacBio::BAM::ExternalResources& DataSet::ExternalResources() const
{ return d_->ExternalResources(); }

inline PacBio::BAM::ExternalResources& DataSet::ExternalResources()
{ return d_->ExternalResources(); }

inline DataSet& DataSet::ExternalResources(const PacBio::BAM::ExternalResources& resources)
{ d_->ExternalResources(resources); return *this; }

inline const PacBio::BAM::Filters& DataSet::Filters() const
{ return d_->Filters(); }

inline PacBio::BAM::Filters& DataSet::Filters()
{ return d_->Filters(); }

inline DataSet& DataSet::Filters(const PacBio::BAM::Filters& filters)
{ d_->Filters(filters); return *this; }

inline const std::string& DataSet::Format() const
{ return d_->Format(); }

inline std::string& DataSet::Format()
{ return d_->Format(); }

inline DataSet& DataSet::Format(const std::string& format)
{ d_->Format(format); return *this; }

inline const PacBio::BAM::DataSetMetadata& DataSet::Metadata() const
{ return d_->Metadata(); }

inline PacBio::BAM::DataSetMetadata& DataSet::Metadata()
{ return d_->Metadata(); }

inline DataSet& DataSet::Metadata(const PacBio::BAM::DataSetMetadata& metadata)
{ d_->Metadata(metadata); return *this; }

inline const std::string& DataSet::MetaType() const
{ return d_->MetaType(); }

inline std::string& DataSet::MetaType()
{ return d_->MetaType(); }

inline DataSet& DataSet::MetaType(const std::string& metatype)
{ d_->MetaType(metatype); return *this; }

inline const std::string& DataSet::ModifiedAt() const
{ return d_->ModifiedAt(); }

inline std::string& DataSet::ModifiedAt()
{ return d_->ModifiedAt(); }

inline DataSet& DataSet::ModifiedAt(const std::string& modifiedAt)
{ d_->ModifiedAt(modifiedAt); return *this; }

inline const std::string& DataSet::Name() const
{ return d_->Name(); }

inline std::string& DataSet::Name()
{ return d_->Name(); }

inline DataSet& DataSet::Name(const std::string& name)
{ d_->Name(name); return *this; }

inline const std::string& DataSet::ResourceId() const
{ return d_->ResourceId(); }

inline std::string& DataSet::ResourceId()
{ return d_->ResourceId(); }

inline DataSet& DataSet::ResourceId(const std::string& resourceId)
{ d_->ResourceId(resourceId); return *this; }

inline const PacBio::BAM::SubDataSets& DataSet::SubDataSets() const
{ return d_->SubDataSets(); }

inline PacBio::BAM::SubDataSets& DataSet::SubDataSets()
{ return d_->SubDataSets(); }

inline DataSet& DataSet::SubDataSets(const PacBio::BAM::SubDataSets& subdatasets)
{ d_->SubDataSets(subdatasets); return *this; }

inline const std::string& DataSet::Tags() const
{ return d_->Tags(); }

inline std::string& DataSet::Tags()
{ return d_->Tags(); }

inline DataSet& DataSet::Tags(const std::string& tags)
{ d_->Tags(tags); return *this; }

inline const std::string& DataSet::TimeStampedName() const
{ return d_->TimeStampedName(); }

inline std::string& DataSet::TimeStampedName()
{ return d_->TimeStampedName(); }

inline DataSet& DataSet::TimeStampedName(const std::string& timeStampedName)
{ d_->TimeStampedName(timeStampedName); return *this; }

inline PacBio::BAM::DataSet::TypeEnum DataSet::Type() const
{ return DataSet::NameToType(TypeName()); }

inline DataSet& DataSet::Type(const DataSet::TypeEnum type)
{ d_->Label(DataSet::TypeToName(type)); return *this; }

inline std::string DataSet::TypeName() const
{ return d_->LocalNameLabel().to_string(); }

inline const std::string& DataSet::UniqueId() const
{ return d_->UniqueId(); }

inline std::string& DataSet::UniqueId()
{ return d_->UniqueId(); }

inline DataSet& DataSet::UniqueId(const std::string& uuid)
{ d_->UniqueId(uuid); return *this; }

inline const std::string& DataSet::Version() const
{ return d_->Version(); }

inline std::string& DataSet::Version()
{ return d_->Version(); }

inline DataSet& DataSet::Version(const std::string& version)
{ d_->Version(version); return *this; }

} // namespace BAM
} // namespace PacBio
