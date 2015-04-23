/* QualityValue.i */
%module PacBioBam
%{
#include <pbbam/QualityValue.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

%ignore PacBio::BAM::QualityValue::operator=;

#ifdef SWIGPYTHON
%rename(__int__) PacBio::BAM::QualityValue::operator uint8_t;
#else // R, C#
%rename(ToInt) PacBio::BAM::QualityValue::operator uint8_t;
#endif

%include <pbbam/QualityValue.h>

