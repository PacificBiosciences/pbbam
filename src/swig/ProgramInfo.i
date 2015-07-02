/* ProgramInfo.i */

%module PacBioBam

%{
#include <pbbam/ProgramInfo.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

%ignore PacBio::BAM::ProgramInfo::ProgramInfo(ProgramInfo&&);
%ignore PacBio::BAM::ProgramInfo::operator=;
%ignore PacBio::BAM::ProgramInfo::ToSam(const ProgramInfo&);    // ignore static method, to allow member

%include <pbbam/ProgramInfo.h>