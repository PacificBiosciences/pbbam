/* SequenceInfo.i */

%module PacBioBam

%{
#include <pbbam/SequenceInfo.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

%ignore PacBio::BAM::SequenceInfo::SequenceInfo(SequenceInfo&&);
%ignore PacBio::BAM::SequenceInfo::operator=;
%ignore PacBio::BAM::SequenceInfo::ToSam(const SequenceInfo&);    // ignore static method, to allow member

%include <pbbam/SequenceInfo.h>