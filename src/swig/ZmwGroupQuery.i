/* ZmwGroupQuery.i */

%module PacBioBam

%{
#include <pbbam/internal/QueryBase.h>
#include <pbbam/ZmwGroupQuery.h>
	
using namespace PacBio;
using namespace PacBio::BAM;
%}

%include <pbbam/internal/QueryBase.h>
%include <pbbam/ZmwGroupQuery.h>
