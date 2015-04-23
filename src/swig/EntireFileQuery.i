/* EntireFileQuery.i */
%module PacBioBam
%{
#include <pbbam/EntireFileQuery.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

HANDLE_STD_EXCEPTION(EntireFileQuery);

%include <pbbam/EntireFileQuery.h>