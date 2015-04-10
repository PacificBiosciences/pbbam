/* BamRecord.i */

%module PacBioBam

%{
#include <pbbam/BamRecord.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

// hide warnings about unused methods
%ignore PacBio::BAM::BamRecord::BamRecord(BamRecordImpl&&); 
%ignore PacBio::BAM::BamRecord::BamRecord(BamRecord&&);
%ignore PacBio::BAM::BamRecord::operator=;

// ignore static methods, to allow member
%ignore PacBio::BAM::BamRecord::Clipped(const BamRecord&, const ClipType, const PacBio::BAM::Position, const PacBio::BAM::Position);
%ignore PacBio::BAM::BamRecord::Mapped(const BamRecord&, const int32_t, const Position, const Strand, const Cigar&, const uint8_t);

%include <pbbam/BamRecord.h>