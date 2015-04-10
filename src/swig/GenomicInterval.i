/* GenomicInterval.i */
%module PacBioBam
%{
#include <pbbam/GenomicInterval.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

%include <pbbam/GenomicInterval.h>