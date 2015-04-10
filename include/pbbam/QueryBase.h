// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

// Author: Derek Barnett

#ifndef QUERYBASE_H
#define QUERYBASE_H

#include "pbbam/BamRecord.h"
#include "pbbam/BamFile.h"

namespace PacBio {
namespace BAM {

class QueryBase;

class QueryIterator
{
public:
    BamRecord& operator*(void);
    BamRecord* operator->(void);
    QueryIterator& operator++(void);
    QueryIterator operator++(int);
    bool operator==(const QueryIterator& other) const;
    bool operator!=(const QueryIterator& other) const;

    QueryIterator(void);
    QueryIterator(QueryBase& parent);

private:
    QueryBase* query_;
    BamRecord record_;
    friend class QueryBase;
};

class QueryConstIterator
{
public:
    const BamRecord& operator*(void) const;
    const BamRecord* operator->(void) const;
    QueryConstIterator& operator++(void);
    QueryConstIterator operator++(int);
    bool operator==(const QueryConstIterator& other) const;
    bool operator!=(const QueryConstIterator& other) const;

    QueryConstIterator(void);
    QueryConstIterator(const QueryBase& parent);

private:
    QueryBase* query_;
    BamRecord record_;
    friend class QueryBase;
};

/// This class provides the base functionality and iterators for querying BAM files.
class PBBAM_EXPORT QueryBase {

public:
    typedef QueryIterator      iterator;
    typedef QueryConstIterator const_iterator;

public:
    /// This enum describes the errors that may be returned by the Error() function.
    enum QueryError
    {
        NoError                   ///< No error occurred.
      , FileOpenError             ///< An error occurred while opening the BAM file.
      , FileMetadataError         ///< An error occurred while reading the BAM file metadata.
      , IndexFileOpenError        ///< An error occurred while opening the index file.
      , IndexFileMetadataError    ///< An error occurred while reading the index file metadata.
      , InitializeQueryError      ///< An error occurred while initializing query (e.g. invalid parameters).
    };


public:
    virtual ~QueryBase(void);

public:

    /// \name Error Handling
    /// \{

    /// \returns the query's error status.
    QueryError Error(void) const;

    /// \returns true if Error() is QueryBase::NoError, else false
    operator bool(void) const;

    /// \}

public:

    /// \name Iterators
    /// \{

    /// \returns an iterator to the beginning of the query results.
    QueryBase::iterator begin(void);

    /// \returns a const_iterator to the beginning of the query results.
    QueryBase::const_iterator begin(void) const;

    /// \returns a const_iterator to the beginning of the query results.
    QueryBase::const_iterator cbegin(void) const;

    /// \returns an iterator marking the end of query results.
    QueryBase::iterator end(void);

    /// \returns a const_iterator marking the end of query results.
    QueryBase::const_iterator end(void) const;

    /// \returns a const_iterator marking the end of query results.
    QueryBase::const_iterator cend(void) const;

    /// \}

protected:
    QueryBase(const BamFile& file);

    /// Primary method for iterating through a query. Derived classes will implement this
    /// method to return
    virtual bool GetNext(BamRecord& x) =0;

protected:
    QueryError error_;
    const BamFile& file_;

    friend class QueryIterator;
    friend class QueryConstIterator;
};

inline QueryBase::iterator QueryBase::begin(void)
{ return QueryBase::iterator(*this); }

inline QueryBase::const_iterator QueryBase::begin(void) const
{ return QueryBase::const_iterator(*this); }

inline QueryBase::const_iterator QueryBase::cbegin(void) const
{ return QueryBase::const_iterator(*this); }

inline QueryBase::iterator QueryBase::end(void)
{ return QueryBase::iterator(); }

inline QueryBase::const_iterator QueryBase::end(void) const
{ return QueryBase::const_iterator(); }

inline QueryBase::const_iterator QueryBase::cend(void) const
{ return QueryBase::const_iterator(); }

inline QueryBase::QueryError QueryBase::Error(void) const
{ return error_; }

inline QueryBase::operator bool(void) const
{ return error_ == QueryBase::NoError; }

// ---------------
// QueryIterator
// ---------------

inline QueryIterator::QueryIterator(void)
    : query_(0)
{ }

inline QueryIterator::QueryIterator(QueryBase& parent)
    : query_(&parent)
    , record_(parent.file_.Header())
{
    if (!(query_->GetNext(record_)))
        query_ = 0;
}

inline BamRecord& QueryIterator::operator*(void)
{ return record_; }

inline BamRecord* QueryIterator::operator->(void)
{ return &(operator*()); }

inline QueryIterator& QueryIterator::operator++(void)
{
    if (!(query_->GetNext(record_)))
        query_ = 0;
    return *this;
}

inline QueryIterator QueryIterator::operator++(int)
{
    QueryIterator result(*this);
    ++(*this);
    return result;
}

inline bool QueryIterator::operator==(const QueryIterator& other) const
{ return query_ == other.query_; }

inline bool QueryIterator::operator!=(const QueryIterator& other) const
{ return !(*this == other); }

// --------------------
// QueryConstIterator
// --------------------

inline const BamRecord& QueryConstIterator::operator*(void) const
{ return record_; }

inline const BamRecord* QueryConstIterator::operator->(void) const
{ return &(operator*()); }

inline QueryConstIterator& QueryConstIterator::operator++(void)
{
    if (!(query_->GetNext(record_)))
        query_ = 0;
    return *this;
}

inline QueryConstIterator QueryConstIterator::operator++(int)
{
    QueryConstIterator result(*this);
    ++(*this);
    return result;
}

inline bool QueryConstIterator::operator==(const QueryConstIterator& other) const
{ return query_ == other.query_; }

inline bool QueryConstIterator::operator!=(const QueryConstIterator& other) const
{ return !(*this == other); }

inline QueryConstIterator::QueryConstIterator(void)
    : query_(0)
{ }

inline QueryConstIterator::QueryConstIterator(const QueryBase& parent)
    : record_(parent.file_.Header())
{
    query_ = const_cast<QueryBase*>(&parent);
    if (!(query_->GetNext(record_)))
        query_ = 0;
}

} // namespace BAM
} // namspace PacBio

#endif // QUERYBASE_H
