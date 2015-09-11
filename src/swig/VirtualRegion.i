/* VirtualRegion.i */

%module PacBioBam

%{
#include <pbbam/virtual/VirtualRegionType.h>
#include <pbbam/virtual/VirtualRegion.h>
#include <map>
#include <vector>
using namespace PacBio;
using namespace PacBio::BAM;
%}

%ignore PacBio::BAM::VirtualRegion::VirtualRegion(VirtualRegion&&);
%ignore PacBio::BAM::VirtualRegion::operator=;

%include <pbbam/virtual/VirtualRegionType.h>
%include <pbbam/virtual/VirtualRegion.h>
