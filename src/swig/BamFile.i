/* BamFile.i */

%module PacBioBam

%{
#include <pbbam/BamFile.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

%ignore PacBio::BAM::BamFile::BamFile(BamFile&&);
%ignore PacBio::BAM::BamFile::operator=;

HANDLE_STD_EXCEPTION(BamFile);
HANDLE_STD_EXCEPTION(EnsurePacBioIndexExists);
HANDLE_STD_EXCEPTION(EnsureStandardIndexExists);
HANDLE_STD_EXCEPTION(ReferenceId);
HANDLE_STD_EXCEPTION(ReferenceLength);
HANDLE_STD_EXCEPTION(ReferenceName);

%include <pbbam/BamFile.h>