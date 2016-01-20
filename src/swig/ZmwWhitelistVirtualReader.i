/* ZmwWhitelistVirtualReader.i */

%module PacBioBam

%{
#include <pbbam/virtual/ZmwWhitelistVirtualReader.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

%include <pbbam/virtual/ZmwWhitelistVirtualReader.h>
