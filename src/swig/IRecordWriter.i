/* IRecordWriter.i */
%module PacBioBam
%{
#include <pbbam/IRecordWriter.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

%include <pbbam/IRecordWriter.h>
