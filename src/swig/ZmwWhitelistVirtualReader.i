/* ZmwWhitelistVirtualReader.i */

%module PacBioBam

%{
#include <pbbam/virtual/ZmwWhitelistVirtualReader.h>
#include <pbbam/virtual/WhitelistedZmwReadStitcher.h>
using namespace PacBio;
using namespace PacBio::BAM;
typedef PacBio::BAM::WhitelistedZmwReadStitcher ZmwWhitelistVirtualReader;
%}

%include <pbbam/virtual/ZmwWhitelistVirtualReader.h>
%include <pbbam/virtual/WhitelistedZmwReadStitcher.h>
typedef PacBio::BAM::WhitelistedZmwReadStitcher ZmwWhitelistVirtualReader;

#ifdef SWIGPYTHON
%pythoncode %{

ZmwWhitelistVirtualReader = WhitelistedZmwReadStitcher

%}
#endif 
