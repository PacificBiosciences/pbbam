/* QueryBase.i */
%module PacBioBam
%{
#include <pbbam/IndexedFastaReader.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}


%include <pbbam/IndexedFastaReader.h>
