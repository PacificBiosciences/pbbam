/* BamRecordImpl.i */

%module PacBioBam

%{
#include <pbbam/BamRecordImpl.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

%ignore PacBio::BAM::BamRecordImpl::BamRecordImpl(BamRecordImpl&&); 
%ignore PacBio::BAM::BamRecordImpl::operator=;

HANDLE_STD_EXCEPTION(CigarData);

%include <pbbam/BamRecordImpl.h>