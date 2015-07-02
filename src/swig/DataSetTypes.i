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
