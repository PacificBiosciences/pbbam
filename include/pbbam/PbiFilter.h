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
// Author: Derek Barnett

#ifndef PBIFILTER_H
#define PBIFILTER_H

#include "pbbam/DataSet.h"
#include "pbbam/PbiBasicTypes.h"
#include "pbbam/PbiIndex.h"
#include <boost/concept_check.hpp>
#include <memory>
#include <string>

namespace PacBio {
namespace BAM {

namespace internal { struct PbiFilterPrivate; }

template<typename T>
struct PbiFilterConcept
{
    BOOST_CONCEPT_USAGE(PbiFilterConcept)
    {
        // PBI filters (to be accepted by PbiIndex & CompositePbiFilter)
        // need only provide this interface:
        //
        //    IndexList Lookup(const PbiIndex& index) const;
        //
        result = filter.Lookup(index);
    }

private:
    T filter;
    PbiIndex index;
    IndexList result;
};

// TODO: clean up this documentation ... move some to a "recipe" section ?

/// PbiIndex filters are designed to be flexible, both built-in and for client-side
/// customization.
///
/// Set Composition - Intersect & Union
///
/// \code
///
/// // (f1 && f2) || f3
/// PbiFilter f1;
/// PbiFilter f2;
/// PbiFilter intersect_f1_f2 = PbiFilter::Intersect(f1, f2);
/// PbiFilter f3;
/// PbiFilter final = PbiFilter::Union(intersect_f1_f2, f3);
///
/// \endcode
///
/// Filter objects used as children of PbiFilter need only provide a method that matches this signature:
///
/// \code {.cpp}
/// IndexList Lookup(const PbiIndex& index) const;
/// \endcode
///
/// This requirement is 'enforced' by using the PbiFilterConcept to require interface without requiring inheritance.
/// This approach allows composition of heterogeneous types without worrying about a class hierarchy, pointer ownership, etc.
///
/// Thus a client application can define a custom filter if the built-in filters do not quite meet their needs:
///
/// \code
/// struct MyCustomFilter
/// {
///     IndexList Lookup(const PbiIndex& index) const
///     {
///         // do index queries on index data
///         // filter, calculate, transform, etc. as needed
///         return IndexList();
///     }
/// };
/// \endcode
///
/// And this filter may be used in PbiFilter composition, or directly to PbiFilterQuery.
///
///
/// \code
/// PbiFilter filter;
/// filter.Add(MyCustomFilter(42));
/// \endcode
///
///
class PBBAM_EXPORT PbiFilter
{

public:
    enum CompositionType
    {
        INTERSECT
      , UNION
    };

public:
    /// \name Set Operations
    /// \{

    /// Constructs a PbiFilter that acts as an intersection of the input filters.
    ///
    /// i.e. records must satisfy all filters
    ///
    /// Exactly equivalent to:
    /// \code{.cpp}
    ///     PbiFilter result{ PbiFilter::INTERSECT };
    ///     result.Add(filters);
    ///     return result;
    /// \endcode
    ///
    static PbiFilter Intersection(const std::vector<PbiFilter>& filters);

    /// Constructs a PbiFilter that acts as an intersection of the input filters.
    ///
    /// i.e. records must satisfy all filters
    ///
    /// Exactly equivalent to:
    /// \code{.cpp}
    ///     PbiFilter result{ PbiFilter::INTERSECT };
    ///     result.Add(std::move(filters));
    ///     return result;
    /// \endcode
    ///
    static PbiFilter Intersection(std::vector<PbiFilter>&& filters);

    /// Constructs a PbiFilter that acts as a union of the input filters.
    ///
    /// i.e. records must satisfy any child filters
    ///
    /// Exactly equivalent to:
    /// \code{.cpp}
    ///     PbiFilter result{ PbiFilter::UNION };
    ///     result.Add(filters);
    ///     return result;
    /// \endcode
    ///
    static PbiFilter Union(const std::vector<PbiFilter>& filters);

    /// Constructs a PbiFilter that acts as a union of the input filters.
    ///
    /// i.e. records must satisfy any child filters
    ///
    /// Exactly equivalent to:
    /// \code{.cpp}
    ///     PbiFilter result{ PbiFilter::UNION };
    ///     result.Add(std::move(filters));
    ///     return result;
    /// \endcode
    ///
    static PbiFilter Union(std::vector<PbiFilter>&& filters);

    /// \}

public:
    /// \name Constructors & Related Methods
    /// \{

    /// Construct a PbiFilter DataSet XML data
    ///
    /// The DataSet XML element contains a Filters element. The resulting PbiFilter from this method
    /// will perform a union of these filters. i.e. records must pass any one of the filters described.
    ///
    ///
    ///
    ///
    static PbiFilter FromDataSet(const DataSet& dataset);

public:

    /// Constructs an empty filter.
    ///
    /// \param[in] type composition type. Any additional child filters added to this composite will be
    ///                 treated according to this type If INTERSECT, a record must match all child filters.
    ///                 If UNION, a record must match any child filter.
    ///
    PbiFilter(const CompositionType type = INTERSECT);

    /// Constructs a composite (of INTERSECT type) with an initial filter.
    ///
    /// \note T must satisfy PbiFilterConcept
    ///
    /// \param[in] filter initial child filter
    ///
    template<typename T>
    PbiFilter(const T& filter);

    /// Constructs a composite (of INTERSECT type) with an initial filter.
    ///
    /// \note T must satisfy PbiFilterConcept
    ///
    /// \param[in] filter initial child filter
    ///
    template<typename T>
    PbiFilter(T&& filter);

    /// Constructs a composite (of INTERSECT type) with a list of initial filters.
    ///
    /// \param[in] filters initial child filters
    ///
    PbiFilter(const std::vector<PbiFilter>& filters);

    /// Constructs a composite (of INTERSECT type) with a list of initial filters.
    ///
    /// \param[in] filters initial child filters
    ///
    PbiFilter(std::vector<PbiFilter>&& filters);

    PbiFilter(const PbiFilter& other);
    PbiFilter(PbiFilter&& other) noexcept;

    PbiFilter& operator=(const PbiFilter& other);
    PbiFilter& operator=(PbiFilter&& other) noexcept;
    ~PbiFilter(void);

    /// \}

public:
    /// \name Composition
    /// \{

    /// Adds a new child filter of type T.
    ///
    /// \note T must satisfy PbiFilterConcept
    ///
    template<typename T>
    PbiFilter& Add(const T& filter);

    /// Adds a new child filter of type T.
    ///
    /// \note T must satisfy PbiFilterConcept
    ///
    template<typename T>
    PbiFilter& Add(T&& filter);

    /// Adds a new child filter.
    ///
    PbiFilter& Add(const PbiFilter& filter);

    /// Adds a new child filter.
    ///
    PbiFilter& Add(PbiFilter&& filter);

    /// Add child filters.
    ///
    PbiFilter& Add(const std::vector<PbiFilter>& filters);

    /// Add child filters.
    ///
    PbiFilter& Add(std::vector<PbiFilter>&& filters);

    /// \}

public:
    /// \name Lookup
    /// \{

    /// Performs the PBI index lookup, combining child results if PbiFilter is a composite.
    ///
    /// \param[in] idx PBI index object
    /// \returns list of \b sorted index positions (record #) that match all filter criteria.
    ///
    IndexList Lookup(const PacBio::BAM::PbiIndex& idx) const;

    /// \}

private:
    std::unique_ptr<internal::PbiFilterPrivate> d_;
};

} // namespace BAM
} // namespace PacBio

#include "pbbam/internal/PbiFilter.inl"
#include "pbbam/PbiFilterTypes.h"

#endif // PBIFILTER_H
