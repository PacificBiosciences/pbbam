/* ReadGroupInfo.i */

%module PacBioBam

%{
#include <pbbam/ReadGroupInfo.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

%ignore PacBio::BAM::ReadGroupInfo::ReadGroupInfo(ReadGroupInfo&&);
%ignore PacBio::BAM::ReadGroupInfo::operator=;
%ignore PacBio::BAM::ReadGroupInfo::ToSam(const ReadGroupInfo&); 

%include <pbbam/ReadGroupInfo.h>