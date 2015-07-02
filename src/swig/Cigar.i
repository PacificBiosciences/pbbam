/* Cigar.i */

%module PacBioBam

%{
#include <pbbam/Cigar.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

%template(CigarOpList) std::vector<PacBio::BAM::CigarOperation>;

%ignore PacBio::BAM::Cigar::Cigar(Cigar&&);
%ignore PacBio::BAM::Cigar::operator=;

%include <pbbam/Cigar.h>
