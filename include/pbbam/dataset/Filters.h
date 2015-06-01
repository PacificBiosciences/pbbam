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

#ifndef FILTERS_H
#define FILTERS_H

#include "pbbam/Config.h"
#include "pbbam/internal/DataSetElement.h"
#include "pbbam/internal/DataSetListElement.h"
#include <string>

namespace PacBio {
namespace BAM {

class PBBAM_EXPORT FilterParameter : public internal::DataSetElement
{
public:
    /// \name Constructors & Related Methods
    /// \{

    FilterParameter(void);

    using DataSetElement::DataSetElement;

    /// \}

public:
    /// \name Attributes

    /// \returns "Name" attribute value
    ///
    const std::string& Name(void) const;

    /// \returns "Value" attribute value
    ///
    const std::string& Value(void) const;

    /// \}

public:
    /// \name Attributes

    /// Sets "Name" attribute.
    ///
    /// \param[in] name
    /// \returns reference to this filter
    ///
    FilterParameter& Name(const std::string& name);

    /// Sets "Value" attribute.
    ///
    /// \param[in] value
    /// \returns reference to this filter
    ///
    FilterParameter& Value(const std::string& value);

    /// \}
};

class PBBAM_EXPORT FilterParameters : public internal::DataSetListElement<FilterParameter>
{
public:
    /// \name Constructors & Related Methods
    /// \{

    FilterParameters(void);

    using DataSetListElement::DataSetListElement;

    /// \}

public:
    /// \name List Modification
    /// \{

    /// Adds \p filter to list
    ///
    /// \param[in] filter
    /// \returns reference to this list
    ///
    void AddParameter(const FilterParameter& param);

    /// Removes \p filter from list
    ///
    /// \param[in] filter
    /// \returns reference to this list
    ///
    void RemoveParameter(const FilterParameter& param);

    /// \}
};

class PBBAM_EXPORT Filter : public internal::DataSetElement
{
public:
    /// \name Constructors & Related Methods
    /// \{

    Filter(void);
    using DataSetElement::DataSetElement;

    /// \}

public:
    /// \name Parameter List Operations
    /// \{

    /// \returns filter parameter list
    const FilterParameters& FilterParameterList(void) const;

    /// \returns number of parameters
    size_t NumFilterParameters(void) const;

    /// \}

public:
    /// \name Parameter List Operations

    /// Adds \p param to filter
    ///
    /// \param[in] param
    /// \returns reference to this filter
    ///
    void AddParameter(const FilterParameter& param);

    /// Removes \p param from filter
    ///
    /// \param[in] param
    /// \returns reference to this filter
    ///
    void RemoveParameter(const FilterParameter& param);

    /// \returns editable parameters list
    FilterParameters& FilterParameterList(void);

    /// \}

};

class PBBAM_EXPORT Filters : public internal::DataSetListElement<Filter>
{
public:
    /// \name Constructors & Related Methods
    /// \{

    Filters(void);

    using DataSetListElement::DataSetListElement;

    /// \}

public:
    /// \name Filter List Operations
    /// \{

    /// Adds \p filter to this list
    ///
    /// \param[in] filter
    /// \returns reference to this filter list
    ///
    void AddFilter(const Filter& filter);

    /// Removes \p filter from this list
    ///
    /// \param[in] param
    /// \returns reference to this filter list
    ///
    void RemoveFilter(const Filter& filter);

    /// \returns number of filters in this list
    size_t NumFilters(void) const;

    /// \}
};

} // namespace BAM
} // namespace PacBio

#endif // FILTERS_H
