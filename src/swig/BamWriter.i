/* BamWriter.i */

%module PacBioBam

%{
#include <pbbam/BamWriter.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}
 
%ignore PacBio::BAM::BamWriter(const BamWriter&);  // copy ctor not used
%ignore PacBio::BAM::BamWriter(BamWriter&&);       // move ctor not used
%ignore PacBio::BAM::BamWriter::operator=;         // assignment operators not used

%include <pbbam/BamWriter.h>
