/* SamTagCodec.i */
%module PacBioBam
%{
#include <pbbam/SamTagCodec.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

%include <pbbam/SamTagCodec.h>