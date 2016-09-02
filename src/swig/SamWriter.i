/* SamWriter.i */

%module PacBioBam

%{
#include <pbbam/SamWriter.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

%ignore PacBio::BAM::SamWriter(const SamWriter&);  // copy ctor not used
%ignore PacBio::BAM::SamWriter(SamWriter&&);       // move ctor not used
%ignore PacBio::BAM::SamWriter::operator=;         // assignment operators not used

%include <pbbam/SamWriter.h>
