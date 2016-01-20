/* IndexedFastaReader.i */

%module PacBioBam

%{
#include <pbbam/IndexedFastaReader.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

%ignore PacBio::BAM::IndexedFastaReader::operator=; // assignment operators not used

%include <pbbam/IndexedFastaReader.h>
