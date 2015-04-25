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


// IEnumerable<BamRecord> interfaces for Queries
%typemap(csinterfaces) PacBio::BAM::QueryBase "global::System.Collections.IEnumerable\n, global::System.Collections.Generic.IEnumerable<PacBio.BAM.BamRecord>\n";
%typemap(cscode) PacBio::BAM::QueryBase %{
    public global::System.Collections.Generic.IEnumerator<PacBio.BAM.BamRecord> GetEnumerator()
    {
        var i = this.cbegin();
        var e = this.cend();
        while (!i.Equals(e))
        {
            yield return i.value();
            i.incr();
        }
    }

    global::System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }
%}


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
