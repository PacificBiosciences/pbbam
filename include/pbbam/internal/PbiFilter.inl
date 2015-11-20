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
//
// File Description
/// \file PbiFilter.inl
/// \brief Inline implementations for the PbiFilter class.
//
// Author: Derek Barnett

#include "pbbam/PbiFilter.h"
#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <vector>

namespace PacBio {
namespace BAM {
namespace internal {

/// \internal
///
/// This class wraps a the basic PBI filter (whether property filter or some operator
/// e.g. union, intersect, etc.). The wrapper allows PbiFilters to hold heterogeneous,
/// recursive filter types - without exposing pointers & worrying about memory ownership
/// issues between client & library.
///
/// Filters can be given by value from client code and we will wrap them for composition.
///
/// \code{.cpp}
///    PbiFilter f1(ZmwFilter(42));
///    PbiFilter f2;
///    f2.Add(QueryLengthFilter(3000, GREATER_THAN_EQUAL));
///    f2.Add(MyApplicationCustomFilter("foo"));
///    PbiFilter intersect = PbiFilter::Intersect(f1, f2);
///    ...
/// \endcode
///
struct FilterWrapper
{
public:
    template<typename T> FilterWrapper(T x);

    FilterWrapper(const FilterWrapper& other);
    FilterWrapper(FilterWrapper&&) noexcept = default;
    FilterWrapper& operator=(const FilterWrapper& other);
    FilterWrapper& operator=(FilterWrapper&&) noexcept = default;
    ~FilterWrapper(void);

public:
    PacBio::BAM::IndexList Lookup(const PacBio::BAM::PbiIndex& idx) const;

private:
    struct WrapperInterface
    {
        virtual ~WrapperInterface(void) = default;
        virtual WrapperInterface* Clone(void) const =0;
        virtual IndexList Lookup(const PacBio::BAM::PbiIndex& idx) const =0;
    };

    template<typename T>
    struct WrapperImpl : public WrapperInterface
    {
        WrapperImpl(T x);
        WrapperImpl(const WrapperImpl& other);
        WrapperInterface* Clone(void) const;
        IndexList Lookup(const PacBio::BAM::PbiIndex& idx) const;
        T data_;
    };

private:
    std::unique_ptr<WrapperInterface> self_;
};

// ---------------
// FilterWrapper
// ---------------

template<typename T>
inline FilterWrapper::FilterWrapper(T x)
    : self_(new WrapperImpl<T>(std::move(x)))
{ }

inline FilterWrapper::FilterWrapper(const FilterWrapper& other)
    : self_(other.self_->Clone())
{ }

inline FilterWrapper& FilterWrapper::operator=(const FilterWrapper& other)
{
    self_.reset(other.self_->Clone());
    return *this;
}

inline FilterWrapper::~FilterWrapper(void) { }

inline IndexList FilterWrapper::Lookup(const PbiIndex& idx) const
{ return self_->Lookup(idx); }

// ----------------
// WrapperImpl<T>
// ----------------

template<typename T>
inline FilterWrapper::WrapperImpl<T>::WrapperImpl(T x)
    : FilterWrapper::WrapperInterface()
    , data_(std::move(x))
{
    BOOST_CONCEPT_ASSERT((PbiFilterConcept<T>));
}

template<typename T>
inline FilterWrapper::WrapperImpl<T>::WrapperImpl(const WrapperImpl& other)
    : FilterWrapper::WrapperInterface()
    , data_(other.data_)
{ }

template<typename T>
inline FilterWrapper::WrapperInterface* FilterWrapper::WrapperImpl<T>::Clone(void) const
{ return new WrapperImpl(*this); }

template<typename T>
inline IndexList FilterWrapper::WrapperImpl<T>::Lookup(const PbiIndex& idx) const
{ return data_.Lookup(idx); }

struct PbiFilterPrivate
{
    PbiFilterPrivate(PbiFilter::CompositionType type)
        : type_(type)
    { }

    template<typename T>
    void Add(T&& filter)
    {
        filters_.emplace_back(std::move(filter));
    }

    std::unique_ptr<internal::PbiFilterPrivate> DeepCopy(void)
    {
        auto copy = std::unique_ptr<PbiFilterPrivate>{ new PbiFilterPrivate{type_} };
        copy->filters_ = this->filters_;
        return copy;
    }

    IndexList Lookup(const PbiIndex& idx) const
    {
        // if no child filters (i.e. this is an empty PbiFilter{ }), return all indices
        if (filters_.empty()) {
            const size_t numRecords = idx.NumReads();
            auto result = IndexList{ };
            result.reserve(numRecords);
            for (size_t i = 0; i < numRecords; ++i)
                result.push_back(i);
            return result;
        }

        // do child filter lookups
        auto resultMap = std::map<size_t, size_t>{ }; // index -> numOccurrences
        for (auto&& filter : filters_) {
            auto filterResult = filter.Lookup(idx);
            for (auto&& i : filterResult)
                ++resultMap[i];
        }

        // extract results, depending on composite type
        // note: if we only have one filter, this will work here too
        //
        const size_t numChildFilters = filters_.size();
        auto result = IndexList{ };
        result.reserve(resultMap.size());
        auto iter = resultMap.cbegin();
        auto end  = resultMap.cend();
        for ( ; iter != end; ++iter ) {
            const size_t i = iter->first;
            const size_t numOccurrences = iter->second;
            if (type_ == PbiFilter::INTERSECT) {
                if (numOccurrences == numChildFilters)
                    result.push_back(i);
            } else if (type_ == PbiFilter::UNION) {
                assert(numOccurrences != 0);
                result.push_back(i);
            }
            else
                assert(false); // invalid composite type
        }

        // result is already sorted and uniq-ed
        return result;
    }

    PbiFilter::CompositionType type_;
    std::vector<FilterWrapper> filters_;
};

} // namespace internal

inline PbiFilter::PbiFilter(const CompositionType type)
    : d_{ new internal::PbiFilterPrivate{ type } }
{ }

template<typename T> inline
PbiFilter::PbiFilter(const T& filter)
    : d_{ new internal::PbiFilterPrivate{ PbiFilter::INTERSECT } }
{
    Add(filter);
}

template<typename T> inline
PbiFilter::PbiFilter(T&& filter)
    : d_{ new internal::PbiFilterPrivate{ PbiFilter::INTERSECT } }
{
    Add(std::move(filter));
}

inline PbiFilter::PbiFilter(const std::vector<PbiFilter>& filters)
    : d_{ new internal::PbiFilterPrivate{ PbiFilter::INTERSECT } }
{
    Add(filters);
}

inline PbiFilter::PbiFilter(std::vector<PbiFilter>&& filters)
    : d_{ new internal::PbiFilterPrivate{ PbiFilter::INTERSECT} }
{
    Add(std::move(filters));
}

inline PbiFilter::PbiFilter(const PbiFilter& other)
    : d_{ other.d_->DeepCopy() }
{ }

inline PbiFilter::PbiFilter(PbiFilter&& other) noexcept
    : d_{ std::move(other.d_) }
{ }

inline PbiFilter& PbiFilter::operator=(const PbiFilter& other)
{
    d_ = other.d_->DeepCopy();
    return *this;
}

inline PbiFilter& PbiFilter::operator=(PbiFilter&& other) noexcept
{
    d_ = std::move(other.d_);
    return *this;
}

inline PbiFilter::~PbiFilter(void) { }

template<typename T>
inline PbiFilter& PbiFilter::Add(const T& filter)
{
    T copy = filter;
    return Add(std::move(copy));
}

template<typename T>
inline PbiFilter& PbiFilter::Add(T&& filter)
{
    d_->Add(std::move(filter));
    return *this;
}

inline PbiFilter& PbiFilter::Add(const PbiFilter& filter)
{
    PbiFilter copy = filter;
    return Add(std::move(copy));
}

inline PbiFilter& PbiFilter::Add(PbiFilter&& filter)
{
    d_->Add(std::move(filter));
    return *this;
}

inline PbiFilter& PbiFilter::Add(const std::vector<PbiFilter>& filters)
{
    std::vector<PbiFilter> copy = filters;
    return Add(std::move(copy));
}

inline PbiFilter& PbiFilter::Add(std::vector<PbiFilter>&& filters)
{
    for (auto&& filter : filters)
        d_->Add(std::move(filter));
    return *this;
}

inline IndexList PbiFilter::Lookup(const PacBio::BAM::PbiIndex& idx) const
{ return d_->Lookup(idx); }

} // namespace BAM
} // namespace PacBio
