// File Description
/// \file PbiFilter.h
/// \brief Defines the PbiFilter class & helper 'concept'.
//
// Author: Derek Barnett

#ifndef PBIFILTER_H
#define PBIFILTER_H

#include <boost/concept_check.hpp>
#include <cstddef>
#include <memory>
#include <string>
#include <tuple>
#include "pbbam/DataSet.h"
#include "pbbam/PbiBasicTypes.h"
#include "pbbam/PbiRawData.h"
#include "pbbam/Unused.h"

namespace PacBio {
namespace BAM {

namespace internal {
struct PbiFilterPrivate;
}

/// \brief The PbiFilterConcept class provides compile-time enforcement of the
///        required interface for PbiFilter's child filters.
///
template <typename T>
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
        INTERSECT,
        UNION
    };

public:
    /// \name Set Operations
    /// \{

    /// \brief Creates a PbiFilter that acts as an intersection of the input
    ///        filters.
    ///
    /// A record must satisfy \b all of this filter's direct "child" filters.
    ///
    /// \param[in] filters  vector of child filters
    /// \returns composite filter
    ///
    static PbiFilter Intersection(std::vector<PbiFilter> filters);

    /// \brief Creates a PbiFilter that acts as a union of the input filters.
    ///
    /// A record must satisfy \b any of this filter's direct "child" filters.
    ///
    /// \param[in] filters  vector of child filters
    /// \returns composite filter
    ///
    static PbiFilter Union(std::vector<PbiFilter> filters);

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
    template <typename T>
    PbiFilter(T filter);

    /// \brief Creates composite filter (of INTERSECT type) with a list of
    ///        initial child filters.
    ///
    /// \param[in] filters initial child filters
    ///
    PbiFilter(std::vector<PbiFilter> filters);

    PbiFilter(const PbiFilter&);
    PbiFilter(PbiFilter&&) noexcept = default;
    PbiFilter& operator=(const PbiFilter&);
    PbiFilter& operator=(PbiFilter&&) noexcept = default;
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
    template <typename T>
    PbiFilter& Add(T filter);

    /// \brief Adds a new child filter.
    ///
    /// \param[in] filter   additional child filter
    /// \returns reference to this filter
    ///
    PbiFilter& Add(PbiFilter filter);

    /// \brief Add child filters.
    ///
    /// \param[in] filters  additional child filters
    /// \returns reference to this filter
    ///
    PbiFilter& Add(std::vector<PbiFilter> filters);

    /// \returns true if this filter has no child filters.
    bool IsEmpty() const;

    /// \returns number of child filters
    size_t NumChildren() const;

    /// \returns filter type (intersect, union)
    CompositionType Type() const;

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

}  // namespace BAM
}  // namespace PacBio

#include "pbbam/PbiFilterTypes.h"
#include "pbbam/internal/PbiFilter.inl"

#endif  // PBIFILTER_H
