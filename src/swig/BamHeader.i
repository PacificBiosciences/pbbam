/* BamHeader.i */

%module PacBioBam

%{
#include <pbbam/BamHeader.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

// Hide warnings about "internal" being a C# reserved word
%warnfilter(314) PacBio::BAM::internal;

%ignore PacBio::BAM::BamHeader::BamHeader(BamHeader&&);      // move ctors not used
%ignore PacBio::BAM::BamHeader::operator=;                   // assignment operators not used

%template(ProgramInfoList)   std::vector<PacBio::BAM::ProgramInfo>;
%template(ReadGroupInfoList) std::vector<PacBio::BAM::ReadGroupInfo>;
%template(SequenceInfoList)  std::vector<PacBio::BAM::SequenceInfo>;

%include <pbbam/BamHeader.h>
