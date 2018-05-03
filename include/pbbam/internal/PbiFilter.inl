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

#include "pbbam/MakeUnique.h"

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
///    PbiFilter f1(PbiZmwFilter(42));
///    PbiFilter f2;
///    f2.Add(PbiQueryLengthFilter(3000, GREATER_THAN_EQUAL));
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
    ~FilterWrapper() = default;

public:
    bool Accepts(const PacBio::BAM::PbiRawData& idx, const size_t row) const;

private:
    struct WrapperInterface
    {
        virtual ~WrapperInterface() = default;
        virtual WrapperInterface* Clone() const =0;
        virtual bool Accepts(const PacBio::BAM::PbiRawData& idx,
                             const size_t row) const =0;
    };

    template<typename T>
    struct WrapperImpl : public WrapperInterface
    {
        WrapperImpl(T x);
        WrapperImpl(const WrapperImpl& other);
        WrapperInterface* Clone() const override;
        bool Accepts(const PacBio::BAM::PbiRawData& idx,
                     const size_t row) const override;
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
    : self_{std::make_unique<WrapperImpl<T>>(std::move(x))}
{ }

inline FilterWrapper::FilterWrapper(const FilterWrapper& other)
    : self_{other.self_->Clone()}
{ }

inline FilterWrapper& FilterWrapper::operator=(const FilterWrapper& other)
{
    self_.reset(other.self_->Clone());
    return *this;
}

inline bool FilterWrapper::Accepts(const PbiRawData& idx, const size_t row) const
{ return self_->Accepts(idx, row); }

// ----------------
// WrapperImpl<T>
// ----------------

template<typename T>
inline FilterWrapper::WrapperImpl<T>::WrapperImpl(T x)
    : FilterWrapper::WrapperInterface{}
    , data_(std::move(x))
{
    BOOST_CONCEPT_ASSERT((PbiFilterConcept<T>));
}

template<typename T>
inline FilterWrapper::WrapperImpl<T>::WrapperImpl(const WrapperImpl& other)
    : FilterWrapper::WrapperInterface{}
    , data_(other.data_)
{ }

template<typename T>
inline FilterWrapper::WrapperInterface* FilterWrapper::WrapperImpl<T>::Clone() const
{ return new WrapperImpl(*this); }

template<typename T>
inline bool FilterWrapper::WrapperImpl<T>::Accepts(const PbiRawData& idx,
                                                   const size_t row) const
{ return data_.Accepts(idx, row); }

struct PbiFilterPrivate
{
    PbiFilterPrivate(PbiFilter::CompositionType type = PbiFilter::INTERSECT)
        : type_{type}
    { }

    template<typename T>
    void Add(T filter)
    {
        filters_.emplace_back(std::move(filter));
    }

    std::unique_ptr<internal::PbiFilterPrivate> DeepCopy()
    {
        auto copy = std::make_unique<PbiFilterPrivate>(type_);
        copy->filters_ = this->filters_;
        return copy;
    }

    bool Accepts(const PbiRawData& idx, const size_t row) const
    {
        // no filter -> accepts every record
        if (filters_.empty())
            return true;

        // intersection of child filters
        if (type_ == PbiFilter::INTERSECT) {
            for (const auto& filter : filters_) {
                if (!filter.Accepts(idx, row))
                    return false; // break early on failure
            }
            return true; // all passed
        }

        // union of child filters
        else if (type_ == PbiFilter::UNION) {
            for (const auto& filter : filters_) {
                if (filter.Accepts(idx, row))
                    return true; // break early on pass
            }
            return false; // none passed
        }

        else
            //assert(false); // invalid composite filter type
            throw std::runtime_error{"invalid composite filter type in PbiFilterPrivate::Accepts"};
    }

    PbiFilter::CompositionType type_;
    std::vector<FilterWrapper> filters_;
};

} // namespace internal

inline PbiFilter::PbiFilter(const CompositionType type)
    : d_{std::make_unique<internal::PbiFilterPrivate>(type) }
{ }

template<typename T> inline
PbiFilter::PbiFilter(T filter)
    : d_{std::make_unique<internal::PbiFilterPrivate>() }
{
    Add(std::move(filter));
}

inline PbiFilter::PbiFilter(std::vector<PbiFilter> filters)
    : d_{std::make_unique<internal::PbiFilterPrivate>() }
{
    Add(std::move(filters));
}

inline PbiFilter::PbiFilter(const PbiFilter& other)
    : d_{ other.d_->DeepCopy() }
{ }

inline PbiFilter& PbiFilter::operator=(const PbiFilter& other)
{
    d_ = other.d_->DeepCopy();
    return *this;
}

inline bool PbiFilter::Accepts(const PacBio::BAM::PbiRawData& idx,
                               const size_t row) const
{ return d_->Accepts(idx, row); }

template<typename T>
inline PbiFilter& PbiFilter::Add(T filter)
{
    d_->Add(std::move(filter));
    return *this;
}

inline PbiFilter& PbiFilter::Add(PbiFilter filter)
{
    d_->Add(std::move(filter));
    return *this;
}

inline PbiFilter& PbiFilter::Add(std::vector<PbiFilter> filters)
{
    for (auto&& filter : filters)
        d_->Add(std::move(filter));
    return *this;
}

inline bool PbiFilter::IsEmpty() const
{ return d_->filters_.empty(); }

inline size_t PbiFilter::NumChildren() const
{ return d_->filters_.size(); }

inline PbiFilter::CompositionType PbiFilter::Type() const
{ return d_->type_; }

} // namespace BAM
} // namespace PacBio
