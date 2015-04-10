/* BamFile.i */
%module PacBioBam
%{
#include <pbbam/BamFile.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

// Python
#ifdef SWIGPYTHON
%rename(__nonzero__) PacBio::BAM::BamFile::operator bool; 
#endif // Python

// R
#ifdef SWIGR
%ignore PacBio::BAM::BamFile::operator bool; // use BamFile$IsOpen() instead
#endif // R

%include <pbbam/BamFile.h>