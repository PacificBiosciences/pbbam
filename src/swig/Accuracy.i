/* Accuracy.i */
%module PacBioBam
%{
#include <pbbam/Accuracy.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

#ifdef SWIGPYTHON
%rename(__int__) PacBio::BAM::Accuracy::operator int;
#else // C#, R
%rename(ToInt) PacBio::BAM::Accuracy::operator int;
#endif

%include <pbbam/Accuracy.h>