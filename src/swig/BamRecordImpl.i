/* BamRecordImpl.i */

%module PacBioBam

%{
#include <pbbam/BamRecordImpl.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

%ignore PacBio::BAM::BamRecordImpl::BamRecordImpl(BamRecordImpl&&); 
%ignore PacBio::BAM::BamRecordImpl::operator=;

%include <pbbam/BamRecordImpl.h>