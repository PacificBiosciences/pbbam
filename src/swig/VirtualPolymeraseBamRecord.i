/* VirtualPolymeraseBamRecord.i */

%module PacBioBam

%{
#include <pbbam/virtual/VirtualRegionType.h>
#include <pbbam/virtual/VirtualRegion.h>
#include <pbbam/virtual/VirtualPolymeraseBamRecord.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

/*%ignore PacBio::BAM::VirtualPolymeraseBamRecord::VirtualPolymeraseBamRecord(const VirtualPolymeraseBamRecord&);*/
%ignore PacBio::BAM::VirtualPolymeraseBamRecord::VirtualPolymeraseBamRecord(VirtualPolymeraseBamRecord&&);
%ignore PacBio::BAM::VirtualPolymeraseBamRecord::operator=;

// disabled - can't get it to work right (at least in Python)
// but the same info is available (& correct) from record.VirtualRegionsTable(regionType)
%ignore PacBio::BAM::VirtualPolymeraseBamRecord::VirtualRegionsMap;

%template(VirtualRegionList) std::vector<PacBio::BAM::VirtualRegion>;
%template(VirtualRegionsMap) std::map<PacBio::BAM::VirtualRegionType, std::vector<PacBio::BAM::VirtualRegion> >;

%include <pbbam/virtual/VirtualPolymeraseBamRecord.h>