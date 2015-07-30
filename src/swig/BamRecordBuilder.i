/* BamRecord.i */

%module PacBioBam

%{
#include <pbbam/BamRecordBuilder.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}


%include <pbbam/BamRecordBuilder.h>
