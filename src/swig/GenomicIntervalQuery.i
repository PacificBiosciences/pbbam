/* GenomicIntervalQuery.i */
%module PacBioBam
%{
#include <pbbam/GenomicIntervalQuery.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

%include <pbbam/GenomicIntervalQuery.h>