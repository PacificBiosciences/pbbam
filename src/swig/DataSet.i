/* DataSet.i */

%module PacBioBam

%{
#include <pbbam/DataSet.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

// move ctors not used
%ignore PacBio::BAM::DataSet::DataSet(DataSet&&);      

// assignment operators not used
%ignore PacBio::BAM::DataSet::operator=;                 

%include <pbbam/DataSet.h>