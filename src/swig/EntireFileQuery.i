/* EntireFileQuery.i */

%module PacBioBam

%{
#include <pbbam/DataSet.h>
#include <pbbam/internal/QueryBase.h>
#include <pbbam/EntireFileQuery.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

%include <pbbam/DataSet.h>
%include <pbbam/internal/QueryBase.h>
%include <pbbam/EntireFileQuery.h>
