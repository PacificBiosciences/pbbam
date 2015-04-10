/* BamWriter.i */
%module PacBioBam
%{
#include <pbbam/BamWriter.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

// Python
#ifdef SWIGPYTHON
%rename(__nonzero__) PacBio::BAM::BamWriter::operator bool; 
#endif // Python

// R
#ifdef SWIGR
%rename(IsOpen) PacBio::BAM::BamWriter::operator bool; 
#endif // R

%include <pbbam/BamWriter.h>