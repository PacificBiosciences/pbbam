/* PbiIndex.i */

%module PacBioBam

%{
#include <pbbam/PbiIndex.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

/*%ignore PacBio::BAM::IndexResultBlock::IndexResultBlock();*/
%ignore PacBio::BAM::IndexResultBlock::IndexResultBlock(size_t, size_t);

%ignore PacBio::BAM::PbiIndex::PbiIndex(PbiIndex&&);      // move ctors not used
%ignore PacBio::BAM::PbiIndex::operator=;                 // assignment operators not used
%ignore PacBio::BAM::PbiIndeX::VirtualFileOffsets; 

%include <pbbam/PbiIndex.h>
