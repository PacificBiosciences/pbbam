/* RecordType.i */

%module PacBioBam

%{
#include <pbbam/RecordType.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

%include <pbbam/RecordType.h>
