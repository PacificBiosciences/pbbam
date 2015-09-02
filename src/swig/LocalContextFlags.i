/* LocalContextFlags.i */

%module PacBioBam

%{
#include <pbbam/LocalContextFlags.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

#ifdef SWIGCSHARP
%ignore operator|(const LocalContextFlags, const LocalContextFlags);
#endif

%include <pbbam/LocalContextFlags.h>
