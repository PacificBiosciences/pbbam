// Copyright (c) 2014, Pacific Biosciences of California, Inc.
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

class PBBAM_EXPORT QueryBase {

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

protected:
    QueryBase(void);
public:
    virtual ~QueryBase(void);

public:

    /// \name Error Handling
    /// \{

    QueryError Error(void) const;
    operator bool(void) const;

    /// \}

    class iterator
    {
    public:

        BamRecord& operator*(void)
        {
            return m_record;
        }

        BamRecord* operator->(void)
        {
            return &(operator*());
        }

        iterator& operator++(void)
        {
            if (!(m_query->GetNext(m_record)))
                m_query = 0;
            return *this;
        }

        iterator operator++(int)
        {
            iterator __t(*this);
            ++(*this);
            return __t;
        }

        inline
        bool operator==(const iterator& other) const
        {
            return m_query == other.m_query;
        }

        inline
        bool operator!=(const iterator& other) const
        {
            return !(*this == other);
        }

    private:
        iterator(void) : m_query(0) { }
        iterator(QueryBase& parent)
            : m_query(&parent)
        {
            if (!(m_query->GetNext(m_record)))
                m_query = 0;
        }

    private:
        QueryBase* m_query;
        BamRecord m_record;
        friend class QueryBase;
    };

    class const_iterator
    {
    public:

        const BamRecord& operator*(void) const
        {
            return m_record;
        }

        const BamRecord* operator->(void) const {
            return &(operator*());
        }

        const_iterator& operator++(void)
        {
            if (!(m_query->GetNext(m_record)))
                m_query = 0;
            return *this;
        }

        const_iterator operator++(int)
        {
            const_iterator __t(*this);
            ++(*this);
            return __t;
        }

        inline
        bool operator==(const const_iterator& other) const
        {
            return m_query == other.m_query;
        }

        inline
        bool operator!=(const const_iterator& other) const
        {
            return !(*this == other);
        }

    private:
        const_iterator(void) : m_query(0) { }
        const_iterator(const QueryBase& parent)
        {
            m_query = const_cast<QueryBase*>(&parent);
            if (!(m_query->GetNext(m_record)))
                m_query = 0;
        }

    private:
        QueryBase* m_query;
        BamRecord m_record;
        friend class QueryBase;
    };

public:

    /// \name Iterators
    /// \{

    QueryBase::iterator begin(void) { return QueryBase::iterator(*this); }
    QueryBase::const_iterator begin(void) const { return QueryBase::const_iterator(*this); }
    QueryBase::const_iterator cbegin(void) const { return QueryBase::const_iterator(*this); }
    QueryBase::iterator end(void) { return QueryBase::iterator(); }
    QueryBase::const_iterator end(void) const { return QueryBase::const_iterator(); }
    QueryBase::const_iterator cend(void) const { return QueryBase::const_iterator(); }

    /// \}

protected:
    virtual bool GetNext(BamRecord& x) =0;

protected:
    QueryError error_;
};

} // namespace BAM
} // namspace PacBio

#endif // QUERYBASE_H
