/* Frames.i */
%module PacBioBam
%{
#include <pbbam/Frames.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

%ignore PacBio::BAM::Frames::Frames(Frames&&);
%ignore PacBio::BAM::Frames::Frames(std::vector<uint16_t>&&);
%ignore PacBio::BAM::Frames::operator=;
%ignore PacBio::BAM::Frames::Data(std::vector<uint16_t>&&);

%template(UInt8List)  std::vector<uint8_t>;
%template(UInt16List) std::vector<uint16_t>;

%include <pbbam/Frames.h>
