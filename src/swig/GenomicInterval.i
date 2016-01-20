/* GenomicInterval.i */

%module PacBioBam

%{
#include <pbbam/GenomicInterval.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

%ignore PacBio::BAM::GenomicInterval::operator=;

%include <pbbam/GenomicInterval.h>
