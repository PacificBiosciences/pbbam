/* VirtualZmwBamRecord.i */

%module PacBioBam

%{
#include <pbbam/virtual/VirtualRegionType.h>
#include <pbbam/virtual/VirtualRegion.h>
#include <pbbam/virtual/VirtualZmwBamRecord.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

%feature("valuewrapper") PacBio::BAM::VirtualZmwBamRecord;

/*%ignore PacBio::BAM::VirtualZmwBamRecord::VirtualZmwBamRecord(const VirtualZmwBamRecord&);*/
%ignore PacBio::BAM::VirtualZmwBamRecord::VirtualZmwBamRecord(VirtualZmwBamRecord&&);
%ignore PacBio::BAM::VirtualZmwBamRecord::operator=;

// disabled - can't get it to work right (at least in Python)
// but the same info is available (& correct) from record.VirtualRegionsTable(regionType)
%ignore PacBio::BAM::VirtualZmwBamRecord::VirtualRegionsMap;

%template(VirtualRegionList) std::vector<PacBio::BAM::VirtualRegion>;
%template(VirtualRegionsMap) std::map<PacBio::BAM::VirtualRegionType, std::vector<PacBio::BAM::VirtualRegion> >;

%include <pbbam/virtual/VirtualZmwBamRecord.h>