/* GroupQueryBase.i */

%module PacBioBam

%{
#include <pbbam/GroupQueryBase.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

%ignore PacBio::BAM::GroupQueryIterator::operator++;
%ignore PacBio::BAM::GroupQueryConstIterator::operator++;

%include <pbbam/GroupQueryBase.h>

%extend PacBio::BAM::GroupQueryIterator
{
    PacBio::BAM::GroupQueryIterator& incr(void)
    { return $self->operator++(); }

    std::vector<PacBio::BAM::BamRecord>* value(void)
    { return $self->operator->(); }
}

%extend PacBio::BAM::GroupQueryConstIterator
{
    PacBio::BAM::GroupQueryConstIterator& incr(void)
    { return $self->operator++(); }

    const std::vector<PacBio::BAM::BamRecord>* value(void) const
    { return $self->operator->(); }
}