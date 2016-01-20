/* Strand.i */

%module PacBioBam

%{
#include <pbbam/Strand.h>	
using namespace PacBio;
using namespace PacBio::BAM;
%}

%include <pbbam/Strand.h>