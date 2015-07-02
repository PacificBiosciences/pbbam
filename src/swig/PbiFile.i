/* PbiFile.i */

%module PacBioBam

%{
#include <pbbam/PbiFile.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

%include <pbbam/PbiFile.h>