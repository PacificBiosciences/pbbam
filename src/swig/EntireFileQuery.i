/* EntireFileQuery.i */
%module PacBioBam
%{

#include <pbbam/EntireFileQuery.h>

using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::staging;
%}

HANDLE_STD_EXCEPTION(EntireFileQuery);

%include <pbbam/EntireFileQuery.h>