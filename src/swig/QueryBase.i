/* QueryBase.i */

%module PacBioBam

%{

#include <pbbam/internal/QueryBase.h>

using namespace PacBio;
using namespace PacBio::BAM;
%}


%ignore PacBio::BAM::QueryIterator::operator++;
%ignore PacBio::BAM::QueryConstIterator::operator++;

%ignore PacBio::BAM::internal::QueryIterator::operator++;
%ignore PacBio::BAM::internal::QueryConstIterator::operator++;

%typemap(csinterfaces) PacBio::BAM::internal::QueryBase<BamRecord>  "global::System.Collections.IEnumerable\n, global::System.Collections.Generic.IEnumerable<PacBio.BAM.BamRecord>\n";
%typemap(cscode) PacBio::BAM::internal::QueryBase<BamRecord>
%{

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

namespace std {
    %template(BamRecordList) std::vector<PacBio::BAM::BamRecord>;
}

%typemap(csinterfaces) PacBio::BAM::internal::QueryBase<std::vector<BamRecord> >  "global::System.Collections.IEnumerable\n, global::System.Collections.Generic.IEnumerable<BamRecordList>\n";
%typemap(cscode) PacBio::BAM::internal::QueryBase<std::vector<BamRecord> >
%{

    public global::System.Collections.Generic.IEnumerator<BamRecordList> GetEnumerator()
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

%include <pbbam/internal/QueryBase.h>

%template(IQuery)      PacBio::BAM::internal::QueryBase<BamRecord>;
%template(IGroupQuery) PacBio::BAM::internal::QueryBase<std::vector<BamRecord> >;

// IEnumerable<BamRecord> interfaces for Queries
%template(BamQueryIteratorBase)       PacBio::BAM::internal::QueryIteratorBase<BamRecord>;
%template(BamGroupQueryIteratorBase)  PacBio::BAM::internal::QueryIteratorBase<std::vector<BamRecord> >;
%template(BamQueryIterator)           PacBio::BAM::internal::QueryIterator<BamRecord>;
%template(BamGroupQueryIterator)      PacBio::BAM::internal::QueryIterator<std::vector<BamRecord> >;
%template(BamQueryConstIterator)      PacBio::BAM::internal::QueryConstIterator<BamRecord>;
%template(BamGroupQueryConstIterator) PacBio::BAM::internal::QueryConstIterator<std::vector<BamRecord> >;

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

%extend PacBio::BAM::internal::QueryIterator<BamRecord>
{
        PacBio::BAM::internal::QueryIterator<BamRecord>& incr(void)
	{ return $self->operator++(); }
	
	PacBio::BAM::BamRecord* value(void)
	{ return $self->operator->(); }
}

%extend PacBio::BAM::internal::QueryConstIterator<BamRecord>
{
        PacBio::BAM::internal::QueryConstIterator<BamRecord>& incr(void)
	{ return $self->operator++(); }
	
	const PacBio::BAM::BamRecord* value(void) const
	{ return $self->operator->(); }
}

%extend PacBio::BAM::internal::QueryIterator<std::vector<BamRecord> >
{
        PacBio::BAM::internal::QueryIterator<std::vector<BamRecord> >& incr(void)
	{ return $self->operator++(); }
	
	std::vector<PacBio::BAM::BamRecord>* value(void)
	{ return $self->operator->(); }
}

%extend PacBio::BAM::internal::QueryConstIterator<std::vector<BamRecord> >
{
        PacBio::BAM::internal::QueryConstIterator<std::vector<BamRecord> >& incr(void)
	{ return $self->operator++(); }
	
	const std::vector<PacBio::BAM::BamRecord>* value(void) const
	{ return $self->operator->(); }
}
