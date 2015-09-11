/* BamFile.i */

%module PacBioBam

%{
#include <pbbam/BamFile.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

#ifdef SWIGR
%ignore PacBio::BAM::BamFile::BamFile(const BamFile&);
#endif 

%ignore PacBio::BAM::BamFile::BamFile(BamFile&&);
%ignore PacBio::BAM::BamFile::operator=;

%include <pbbam/BamFile.h>
