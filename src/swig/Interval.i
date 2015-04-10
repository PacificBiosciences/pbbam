/* Interval.i */
%module PacBioBam
%{
#include <pbbam/Interval.h>
#include <pbbam/Position.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

%include <pbbam/Interval.h>

%template(PositionInterval) PacBio::BAM::Interval<PacBio::BAM::Position>;
