/* LocalContextFlags.i */

%module PacBioBam

%{
#include <pbbam/LocalContextFlags.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

%include <pbbam/LocalContextFlags.h>
