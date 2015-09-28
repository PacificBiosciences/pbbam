/* Accuracy.i */

%module PacBioBam

%{
#include <pbbam/Accuracy.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

#ifdef SWIGPYTHON
%rename(__float__) PacBio::BAM::Accuracy::operator float;
#else // C#, R
%rename(ToFloat) PacBio::BAM::Accuracy::operator float;
#endif

%include <pbbam/Accuracy.h>