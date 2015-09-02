/* DataSetTypes.i */

%module PacBioBam

%{
#include <pbbam/internal/DataSetElement.h>
#include <pbbam/internal/DataSetListElement.h>
#include <pbbam/internal/DataSetBaseTypes.h>
#include <pbbam/DataSetTypes.h>

using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
%}

%ignore PacBio::BAM::internal::DataSetElement::DataSetElement(DataSetElement&&);
%ignore PacBio::BAM::internal::DataSetElement::operator=;
%ignore PacBio::BAM::internal::DataSetElement::operator[];
/*%rename(__getitem__) PacBio::BAM::internal::DataSetElement::operator[];*/

%ignore PacBio::BAM::internal::XmlName::XmlName(XmlName&&);
%ignore PacBio::BAM::internal::XmlName::operator=;

#ifdef SWIGCSHARP

// ignore non-const accessors
%ignore PacBio::BAM::DataSetBase::ExternalResources();
%ignore PacBio::BAM::DataSetBase::Filters();
%ignore PacBio::BAM::DataSetBase::Metadata();
%ignore PacBio::BAM::DataSetBase::Namespaces();
%ignore PacBio::BAM::DataSetBase::SubDataSets();
%ignore PacBio::BAM::DataSetMetadata::NumRecords();
%ignore PacBio::BAM::DataSetMetadata::Provenance();
%ignore PacBio::BAM::DataSetMetadata::TotalLength();
%ignore PacBio::BAM::ExternalResource::ExternalResources();
%ignore PacBio::BAM::Filter::Properties();
%ignore PacBio::BAM::Property::Name();
%ignore PacBio::BAM::Property::Operator();
%ignore PacBio::BAM::Property::Value();
%ignore PacBio::BAM::Provenance::CreatedBy();
%ignore PacBio::BAM::Provenance::CommonServicesInstanceId();
%ignore PacBio::BAM::Provenance::CreatorUserId();
%ignore PacBio::BAM::Provenance::ParentJobId();
%ignore PacBio::BAM::Provenance::ParentTool();
%ignore PacBio::BAM::internal::BaseEntityType::Description();
%ignore PacBio::BAM::internal::BaseEntityType::Extensions();
%ignore PacBio::BAM::internal::BaseEntityType::Format();
%ignore PacBio::BAM::internal::BaseEntityType::ModifiedAt();
%ignore PacBio::BAM::internal::BaseEntityType::Name();
%ignore PacBio::BAM::internal::BaseEntityType::ResourceId();
%ignore PacBio::BAM::internal::BaseEntityType::Tags();
%ignore PacBio::BAM::internal::BaseEntityType::Version();
%ignore PacBio::BAM::internal::DataEntityType::Checksum();
%ignore PacBio::BAM::internal::DataEntityType::EncodedValue();
%ignore PacBio::BAM::internal::DataEntityType::MetaType();
%ignore PacBio::BAM::internal::DataEntityType::SimpleValue();
%ignore PacBio::BAM::internal::DataEntityType::TimeStampedName();
%ignore PacBio::BAM::internal::DataEntityType::UniqueId();
%ignore PacBio::BAM::internal::DataEntityType::ValueDataType();
%ignore PacBio::BAM::internal::DataSetElement::Attribute(const std::string&);
%ignore PacBio::BAM::internal::DataSetElement::Attributes();
%ignore PacBio::BAM::internal::DataSetElement::Children();
%ignore PacBio::BAM::internal::DataSetElement::ChildText(const std::string&);
%ignore PacBio::BAM::internal::DataSetElement::CreatedAt();
%ignore PacBio::BAM::internal::DataSetElement::Text();
%ignore PacBio::BAM::internal::IndexedDataType::FileIndices();
%ignore PacBio::BAM::internal::StrictEntityType::MetaType();
%ignore PacBio::BAM::internal::StrictEntityType::TimeStampedName();
%ignore PacBio::BAM::internal::StrictEntityType::UniqueId();

// disable operator(s)
%ignore PacBio::BAM::DataSetMetadata::operator+=;
%ignore PacBio::BAM::ExternalResources::operator+=;
%ignore PacBio::BAM::Filters::operator+=;
%ignore PacBio::BAM::DataSetBase::operator+=;
%ignore PacBio::BAM::SubDataSets::operator+=;

#endif // C#

%include <pbbam/internal/DataSetElement.h>

%ignore PacBio::BAM::internal::DataSetElementList::operator[];
%ignore PacBio::BAM::internal::DataSetListIterator::operator++;
%ignore PacBio::BAM::internal::DataSetListConstIterator::operator++;

%include <pbbam/internal/DataSetListElement.h>

%template(ExtensionListElement)        PacBio::BAM::internal::DataSetListElement<PacBio::BAM::ExtensionElement>;
%template(ExternalResourceListElement) PacBio::BAM::internal::DataSetListElement<PacBio::BAM::ExternalResource>;
%template(FileIndexListElement)        PacBio::BAM::internal::DataSetListElement<PacBio::BAM::FileIndex>;
%template(FilterListElement)           PacBio::BAM::internal::DataSetListElement<PacBio::BAM::Filter>;
%template(PropertyListElement)         PacBio::BAM::internal::DataSetListElement<PacBio::BAM::Property>;
%template(SubDataSetListElement)       PacBio::BAM::internal::DataSetListElement<PacBio::BAM::DataSetBase>;

%extend PacBio::BAM::internal::DataSetListElement<PacBio::BAM::ExtensionElement> {
    PacBio::BAM::ExtensionElement& __getitem__(unsigned int i)       { return $self->Child<PacBio::BAM::ExtensionElement>(i); }
    PacBio::BAM::ExtensionElement& __getitem__(const std::string& s) { return $self->Child<PacBio::BAM::ExtensionElement>(s); }
}
%extend PacBio::BAM::internal::DataSetListElement<PacBio::BAM::ExternalResource> {
    PacBio::BAM::ExternalResource& __getitem__(unsigned int i)       { return $self->Child<PacBio::BAM::ExternalResource>(i); }
    PacBio::BAM::ExternalResource& __getitem__(const std::string& s) { return $self->Child<PacBio::BAM::ExternalResource>(s); }
}
%extend PacBio::BAM::internal::DataSetListElement<PacBio::BAM::FileIndex> {
    PacBio::BAM::FileIndex& __getitem__(unsigned int i)       { return $self->Child<PacBio::BAM::FileIndex>(i);}
    PacBio::BAM::FileIndex& __getitem__(const std::string& s) { return $self->Child<PacBio::BAM::FileIndex>(s);}
}
%extend PacBio::BAM::internal::DataSetListElement<PacBio::BAM::Filter> {
    PacBio::BAM::Filter& __getitem__(unsigned int i)       { return $self->Child<PacBio::BAM::Filter>(i); }
    PacBio::BAM::Filter& __getitem__(const std::string& s) { return $self->Child<PacBio::BAM::Filter>(s); }
}
%extend PacBio::BAM::internal::DataSetListElement<PacBio::BAM::Property> {
    PacBio::BAM::Property& __getitem__(unsigned int i)       { return $self->Child<PacBio::BAM::Property>(i); }
    PacBio::BAM::Property& __getitem__(const std::string& s) { return $self->Child<PacBio::BAM::Property>(s); }
}
%extend PacBio::BAM::internal::DataSetListElement<PacBio::BAM::DataSetBase> {
    PacBio::BAM::DataSetBase& __getitem__(unsigned int i)       { return $self->Child<PacBio::BAM::DataSetBase>(i); }
    PacBio::BAM::DataSetBase& __getitem__(const std::string& s) { return $self->Child<PacBio::BAM::DataSetBase>(s); }
}

%include <pbbam/internal/DataSetBaseTypes.h>
%include <pbbam/DataSetTypes.h>
