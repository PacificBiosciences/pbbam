/* BamRecord.i */

%module PacBioBam

%{
#include <pbbam/BamRecordBuilder.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

%ignore PacBio::BAM::BamRecordBuilder::BamRecordBuilder(BamRecordBuilder&&);      // move ctors not used
%ignore PacBio::BAM::BamRecordBuilder::operator=;

%ignore PacBio::BAM::BamRecordBuilder::Reset(BamRecord&&);
%ignore PacBio::BAM::BamRecordBuilder::Cigar(PacBio::BAM::Cigar&&);
%ignore PacBio::BAM::BamRecordBuilder::Tags(TagCollection&&);

%include <pbbam/BamRecordBuilder.h>
