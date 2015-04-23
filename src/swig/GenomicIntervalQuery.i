/* GenomicIntervalQuery.i */
%module PacBioBam
%{
#include <pbbam/GenomicIntervalQuery.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

HANDLE_STD_EXCEPTION(GenomicIntervalQuery);
HANDLE_STD_EXCEPTION(Interval);

%include <pbbam/GenomicIntervalQuery.h>