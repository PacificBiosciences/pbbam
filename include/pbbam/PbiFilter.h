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
/// \file PbiFilter.h
/// \brief Defines the PbiFilter class & helper 'concept'.
//
// Author: Derek Barnett

#ifndef PBIFILTER_H
#define PBIFILTER_H

#include "pbbam/DataSet.h"
#include "pbbam/PbiBasicTypes.h"
#include "pbbam/PbiIndex.h"
#include "pbbam/PbiRawData.h"
#include "pbbam/Unused.h"
#include <boost/concept_check.hpp>
#include <cstddef>
#include <memory>
#include <string>
#include <tuple>

namespace PacBio {
namespace BAM {

namespace internal { struct PbiFilterPrivate; }

/// \brief The PbiFilterConcept class provides compile-time enforcement of the
///        required interface for PbiFilter's child filters.
///
template<typename T>
struct PbiFilterConcept
{
    BOOST_CONCEPT_USAGE(PbiFilterConcept)
    {
        // All PBI filters (built-in or client-define) need only provide this
        // interface:
        //
        //    bool Accepts(const PbiRawData& index, const size_t row) const;
        //
        PbiRawData index;
        auto result = filter.Accepts(index, 0);
        UNUSED(result);
    }

private:
    T filter;
//    PbiRawData index;
};

/// \brief The PbiFilter class provides a mechanism for performing PBI-enabled
///        lookups.
///
/// The PbiFilter API is designed to be flexible, both built-in and for
/// client-side customization. Built-in filters are provided for common queries,
/// and client code can define and use custom filters as well. More complex
/// filtering rules can be constructed via composition of simpler child filters.
///
/// Filter objects used as children of PbiFilter need only provide a method that
/// matches this signature:
///
/// \include code/PbiFilter_Interface.txt
///
/// This requirement is enforced internally, using the PbiFilterConcept to
/// require a compatible interface without requiring inheritance. This approach
/// allows composition of heterogeneous filter types without worrying about a
/// class hierarchy, pointer ownership across library/client boundaries, etc.
///
/// Thus a client application can define a custom filter if the built-in filters
/// do not quite meet requirements. This filter may then be used in further
/// PbiFilter composition, or directly to PbiFilterQuery
///
/// \include code/PbiFilter_CustomFilter.txt
///
/// As mentioned above, complex filters can be built up using multiple "child"
/// filters. These complex filters are constructed by using either
/// PbiFilter::Union (logical-OR over all direct children) or
/// PbiFilter::Intersection (logical-AND over direct children).
///
/// \include code/PbiFilter_Composition.txt
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

    /// \brief Creates a PbiFilter that acts as intersection of the input
    ///        filters.
    ///
    /// A record must satisfy \b all of this filter's direct "child" filters.
    ///
    /// Equivalent to:
    /// \include code/PbiFilter_Intersection_Copy.txt
    ///
    /// \param[in] filters  vector of child filters
    /// \returns composite filter
    ///
    static PbiFilter Intersection(const std::vector<PbiFilter>& filters);

    /// \brief Creates a PbiFilter that acts as an intersection of the input
    ///        filters.
    ///
    /// A record must satisfy \b all of this filter's direct "child" filters.
    ///
    /// Equivalent to:
    /// \include code/PbiFilter_Intersection_Move.txt
    ///
    /// \param[in] filters  vector of child filters
    /// \returns composite filter
    ///
    static PbiFilter Intersection(std::vector<PbiFilter>&& filters);

    /// \brief Creates a PbiFilter that acts as a union of the input filters.
    ///
    /// A record must satisfy \b any of this filter's direct "child" filters.
    ///
    /// Equivalent to:
    /// \include code/PbiFilter_Union_Copy.txt
    ///
    /// \param[in] filters  vector of child filters
    /// \returns composite filter
    ///
    static PbiFilter Union(const std::vector<PbiFilter>& filters);

    /// \brief Creates a PbiFilter that acts as a union of the input filters.
    ///
    /// A record must satisfy \b any of this filter's direct "child" filters.
    ///
    /// Equivalent to:
    /// \include code/PbiFilter_Union_Move.txt
    ///
    /// \param[in] filters  vector of child filters
    /// \returns composite filter
    ///
    static PbiFilter Union(std::vector<PbiFilter>&& filters);

    /// \}

public:
    /// \name Constructors & Related Methods
    /// \{

    /// \brief Creates a PbiFilter from a %DataSet's described filters.
    ///
    /// A DataSet may contain a Filters element, itself a list of Filter
    /// elements. Each Filter element will contain a Properties element, itself
    /// a list of Property elements.
    ///
    /// The Filters hierarchy looks like this (in its XML output):
    /// \verbinclude examples/plaintext/PbiFilter_DataSetXmlFilters.txt
    ///
    /// The resulting PbiFilter represents a union over all Filter elements,
    /// with each Filter element requiring an intersection of all of its
    /// Property criteria. These Property elements are mapped to built-in PBI
    /// filter types. To use the labels in the example XML above, the filter
    /// created here is equivalent to:
    ///
    /// (A && B) || (C && D)
    ///
    /// If a DataSet lacks any Filters, then an empty PbiFilter will be created
    /// - corresponding to the dataset's entire contents.
    ///
    /// \param[in] dataset  maybe containing filters
    /// \returns composite filter
    ///
    static PbiFilter FromDataSet(const DataSet& dataset);

public:

    /// \brief Creates an empty filter.
    ///
    /// \note An empty filter will result in all records being returned, e.g.
    ///       for query iteration.
    ///
    /// \param[in] type composition type. Any additional child filters added to
    ///                 this composite will be treated according to this type.
    ///                 If INTERSECT, a record must match all child filters. If
    ///                 UNION, a record must match any child filter.
    ///
    PbiFilter(const CompositionType type = INTERSECT);

    /// \brief Creates a composite filter (of INTERSECT type) with an initial
    ///        child filter.
    ///
    /// \note T must satisfy PbiFilterConcept
    ///
    /// \param[in] filter initial child filter
    ///
    template<typename T>
    PbiFilter(const T& filter);

    /// \brief Creates a composite filter (of INTERSECT type) with an initial
    ///        child filter.
    ///
    /// \note T must satisfy PbiFilterConcept
    ///
    /// \param[in] filter initial child filter
    ///
    template<typename T>
    PbiFilter(T&& filter);

    /// \brief Creates a composite filter (of INTERSECT type) with a list of
    ///        initial child filters.
    ///
    /// \param[in] filters initial child filters
    ///
    PbiFilter(const std::vector<PbiFilter>& filters);

    /// \brief Creates composite filter (of INTERSECT type) with a list of
    ///        initial child filters.
    ///
    /// \param[in] filters initial child filters
    ///
    PbiFilter(std::vector<PbiFilter>&& filters);

    PbiFilter(const PbiFilter& other);
    PbiFilter(PbiFilter&& other) noexcept = default;
    PbiFilter& operator=(const PbiFilter& other);
    PbiFilter& operator=(PbiFilter&& other) noexcept = default;
    ~PbiFilter() = default;

    /// \}

public:
    /// \name Composition
    /// \{

    /// \brief Adds a new child filter of type T.
    ///
    /// \param[in] filter   additional child filter. Type T must satisfy
    ///                     PbiFilterConcept.
    /// \returns reference to this filter
    ///
    template<typename T>
    PbiFilter& Add(const T& filter);

    /// \brief Adds a new child filter of type T.
    ///
    /// \param[in] filter   additional child filter. Type T must satisfy
    ///                     PbiFilterConcept.
    /// \returns reference to this filter
    ///
    template<typename T>
    PbiFilter& Add(T&& filter);

    /// \brief Adds a new child filter.
    ///
    /// \param[in] filter   additional child filter
    /// \returns reference to this filter
    ///
    PbiFilter& Add(const PbiFilter& filter);

    /// \brief Adds a new child filter.
    ///
    /// \param[in] filter   additional child filter
    /// \returns reference to this filter
    ///
    PbiFilter& Add(PbiFilter&& filter);

    /// \brief Add child filters.
    ///
    /// \param[in] filters  additional child filters
    /// \returns reference to this filter
    ///
    PbiFilter& Add(const std::vector<PbiFilter>& filters);

    /// \brief Add child filters.
    ///
    /// \param[in] filters  additional child filters
    /// \returns reference to this filter
    ///
    PbiFilter& Add(std::vector<PbiFilter>&& filters);

    /// \returns true if this filter has no child filters.
    bool IsEmpty() const;

    /// \}

public:
    /// \name Lookup
    /// \{

    /// \brief Performs the PBI index lookup, combining child results a
    ///        composite filter.
    ///
    /// \param[in] idx  PBI (raw) index object
    /// \param[in] row  record number in %BAM/PBI files
    ///
    /// \returns true if record at \p row passes this filter criteria, 
    ///          including children (if any)
    ///
    bool Accepts(const BAM::PbiRawData& idx, const size_t row) const;

    /// \}

private:
    std::unique_ptr<internal::PbiFilterPrivate> d_;
};

} // namespace BAM
} // namespace PacBio

#include "pbbam/internal/PbiFilter.inl"
#include "pbbam/PbiFilterTypes.h"

#endif // PBIFILTER_H
