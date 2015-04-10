/* BamTagCodec.i */
%module PacBioBam
%{
#include <pbbam/BamTagCodec.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

%include <pbbam/BamTagCodec.h>