/* GroupQuery.i */

%module PacBioBam

%{
#include <pbbam/GroupQuery.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

%include <pbbam/GroupQuery.h>