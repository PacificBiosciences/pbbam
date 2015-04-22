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

// Author: Yuan Li
// TODO: Up to Derek's decision. This class mostly references 
// QueryBase.h. We may make QueryBase a template class and make
// GroupQueryBase a specialization of the template.

#ifndef _GROUP_QUERY_BASE_H_
#define _GROUP_QUERY_BASE_H_
#include "pbbam/QueryBase.h"
#include "pbbam/BamRecord.h"
#include <memory>
#include <vector>

namespace PacBio {
namespace BAM {

class GroupQueryBase;

class GroupQueryIterator
{
public:
    std::vector<BamRecord> & operator* (void);
    std::vector<BamRecord> * operator-> (void);
    GroupQueryIterator& operator++ (void);
    GroupQueryIterator  operator++ (int);

    bool operator== (const GroupQueryIterator & other) const;
    bool operator!= (const GroupQueryIterator & other) const;

    GroupQueryIterator(void);
    GroupQueryIterator(GroupQueryBase & parent);

private:
    GroupQueryBase * query_;
    std::vector<BamRecord> records_;
    friend class GroupQueryBase;
};

class GroupQueryConstIterator
{
public:
    const std::vector<BamRecord>& operator*(void) const;
    const std::vector<BamRecord>* operator->(void) const;
    GroupQueryConstIterator& operator++(void);
    GroupQueryConstIterator operator++(int);
    bool operator==(const GroupQueryConstIterator& other) const;
    bool operator!=(const GroupQueryConstIterator& other) const;

    GroupQueryConstIterator(void);
    GroupQueryConstIterator(const GroupQueryBase& parent);

private:
    GroupQueryBase* query_;
    std::vector<BamRecord> records_;
    friend class GroupQueryBase;
};


class PBBAM_EXPORT GroupQueryBase
{
public:
    typedef GroupQueryIterator iterator;

protected:
    QueryBase::QueryError error_;
    BamFile file_;

public:
    QueryBase::QueryError Error(void) const;
    operator bool(void) const;
    virtual ~GroupQueryBase(void);

public:
    GroupQueryBase::iterator begin(void);
    GroupQueryBase::iterator end(void);

protected:
    GroupQueryBase(const BamFile & file);
    virtual bool GetNext(std::vector<BamRecord>& records) = 0;

    friend class GroupQueryIterator;
    friend class GroupQueryConstIterator;
};

inline GroupQueryBase::iterator GroupQueryBase::begin(void)
{ return GroupQueryBase::iterator(*this); }

inline GroupQueryBase::iterator GroupQueryBase::end(void)
{ return GroupQueryBase::iterator(); }

inline QueryBase::QueryError GroupQueryBase::Error(void) const
{ return error_;}

inline GroupQueryBase::operator bool(void) const
{ return error_ == QueryBase::NoError; }

inline GroupQueryBase::GroupQueryBase(const BamFile & file)
    : error_(QueryBase::NoError)
    , file_(file)
{ }

inline GroupQueryBase::~GroupQueryBase(void) { }

// -------------------
// GroupQueryIterator
// -------------------

inline GroupQueryIterator::GroupQueryIterator(void): query_(0) {}

inline GroupQueryIterator::GroupQueryIterator(GroupQueryBase & parent)
    : query_(& parent)
    , records_()
{
    if (!(query_->GetNext(records_)))
        query_ = 0;
}

inline std::vector<BamRecord>& GroupQueryIterator::operator* (void)
{ return records_; }

inline std::vector<BamRecord>* GroupQueryIterator::operator-> (void)
{ return &(operator*()); }

inline GroupQueryIterator& GroupQueryIterator::operator++ (void)
{
    if (!(query_->GetNext(records_)))
        query_ = 0;
    return *this;
}

inline GroupQueryIterator GroupQueryIterator::operator++ (int)
{
    GroupQueryIterator result(*this);
    ++(*this);
    return result;
}

inline bool GroupQueryIterator::operator==(const GroupQueryIterator& other) const
{ return query_ == other.query_; }

inline bool GroupQueryIterator::operator!=(const GroupQueryIterator& other) const
{ return !(*this == other); }


// -------------------
// GroupQueryConstIterator
// -------------------

inline const std::vector<BamRecord>& GroupQueryConstIterator::operator*(void) const
{ return records_; }

inline const std::vector<BamRecord>* GroupQueryConstIterator::operator->(void) const
{ return &(operator*()); }

inline GroupQueryConstIterator& GroupQueryConstIterator::operator++(void)
{
    if (!(query_->GetNext(records_)))
        query_ = 0;
    return *this;
}

inline GroupQueryConstIterator GroupQueryConstIterator::operator++(int)
{
    GroupQueryConstIterator result(*this);
    ++(*this);
    return result;
}

inline bool GroupQueryConstIterator::operator==(const GroupQueryConstIterator& other) const
{ return query_ == other.query_; }

inline bool GroupQueryConstIterator::operator!=(const GroupQueryConstIterator& other) const
{ return !(*this == other); }

inline GroupQueryConstIterator::GroupQueryConstIterator(void): query_(0) { }

inline GroupQueryConstIterator::GroupQueryConstIterator(const GroupQueryBase& parent)
    : query_(const_cast<GroupQueryBase*>(&parent))
    , records_()
{
    if (!(query_->GetNext(records_)))
        query_ = 0;
}

} // namespace BAM
} // namespace PacBio

#endif // _GROUP_QUERY_BASE_H_
