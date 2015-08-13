/* BamRecord.i */

%module PacBioBam

%{
#include <pbbam/BamRecord.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

// Hide warnings about "internal" being a C# reserved word
%warnfilter(314) PacBio::BAM::internal;

// hide warnings about unused methods
%ignore PacBio::BAM::BamRecord::BamRecord(BamRecordImpl&&);
%ignore PacBio::BAM::BamRecord::BamRecord(BamRecord&&);
%ignore PacBio::BAM::BamRecord::operator=;

// ignore static methods, to allow member
%ignore PacBio::BAM::BamRecord::Clipped(const BamRecord&, const ClipType, const PacBio::BAM::Position, const PacBio::BAM::Position);
%ignore PacBio::BAM::BamRecord::Mapped(const BamRecord&, const int32_t, const Position, const Strand, const Cigar&, const uint8_t);

// C# gets confused by the const and nonconst overloads
%ignore PacBio::BAM::BamRecord::Impl() const;

#ifdef SWIGR
%rename("EncodedPkmean") PacBio::BAM::BamRecord::Pkmean(const std::vector<uint16_t>&);
%rename("EncodedPkmid")  PacBio::BAM::BamRecord::Pkmid(const std::vector<uint16_t>&);
#endif // SWIGR

HANDLE_STD_EXCEPTION(CigarData);

%include <pbbam/BamRecord.h>
