/* DataSet.i */

%module PacBioBam

%{
#include <pbbam/DataSet.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

// move ctors not used
%ignore PacBio::BAM::DataSet::DataSet(DataSet&&);      

// assignment operators not used
%ignore PacBio::BAM::DataSet::operator=;                 

#ifdef SWIGCSHARP

// ignore non-const accessors
%ignore PacBio::BAM::DataSet::Attribute(const std::string&);
%ignore PacBio::BAM::DataSet::CreatedAt();
%ignore PacBio::BAM::DataSet::Extensions();
%ignore PacBio::BAM::DataSet::ExternalResources();
%ignore PacBio::BAM::DataSet::Filters();
%ignore PacBio::BAM::DataSet::Format();
%ignore PacBio::BAM::DataSet::Metadata();
%ignore PacBio::BAM::DataSet::MetaType();
%ignore PacBio::BAM::DataSet::ModifiedAt();
%ignore PacBio::BAM::DataSet::Name();
%ignore PacBio::BAM::DataSet::Namespaces();
%ignore PacBio::BAM::DataSet::ResourceId();
%ignore PacBio::BAM::DataSet::SubDataSets();
%ignore PacBio::BAM::DataSet::Tags();
%ignore PacBio::BAM::DataSet::TimeStampedName();
%ignore PacBio::BAM::DataSet::UniqueId();
%ignore PacBio::BAM::DataSet::Version();

// disable operator(s)
%ignore PacBio::BAM::DataSet::operator+=;

#endif // C#

%include <pbbam/DataSet.h>
