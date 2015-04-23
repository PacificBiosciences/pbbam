/* QueryBase.i */
%module PacBioBam
%{
#include <pbbam/QueryBase.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

%ignore PacBio::BAM::QueryIterator::operator++;
%ignore PacBio::BAM::QueryConstIterator::operator++;

// Python
#ifdef SWIGPYTHON
%rename(__nonzero__) PacBio::BAM::QueryBase::operator bool; 
#elif defined(SWIGR)
%rename(isTRUE) PacBio::BAM::QueryBase::operator bool;
#else // C#
%rename(IsTrue) PacBio::BAM::QueryBase::operator bool;
#endif

%include <pbbam/QueryBase.h>

// Iterator API
#ifdef SWIGPYTHON
%pythoncode %{
def Iterate(c):
	i = c.begin()
	e = c.end()
	while i != e:
		yield i.value()
		i.incr()
%}
#endif

%extend PacBio::BAM::QueryIterator
{		
	PacBio::BAM::QueryIterator& incr(void) 
	{ return $self->operator++(); }
	
	PacBio::BAM::BamRecord* value(void) 
	{ return $self->operator->(); }
}

%extend PacBio::BAM::QueryConstIterator
{	
	PacBio::BAM::QueryConstIterator& incr(void) 
	{ return $self->operator++(); }
	
	const PacBio::BAM::BamRecord* value(void) const 
	{ return $self->operator->(); }
}