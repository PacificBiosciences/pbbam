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

#include "pbbam/dataset/Filters.h"
using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

// shared empty, "null" lists
static
const FilterParameters& NullParams(void)
{
    static const FilterParameters empty;
    return empty;
}

} // namespace internal
} // namespace BAM
} // namespace PacBio

// --------------------------------
// FilterParameter implementation
// --------------------------------

FilterParameter::FilterParameter(void)
    : DataSetElement("Parameter")
{ }

const string &FilterParameter::Name(void) const
{ return Attribute("Name"); }

FilterParameter& FilterParameter::Name(const string& name)
{ Attribute("Name", name); return *this; }

const string& FilterParameter::Value(void) const
{ return Attribute("Value"); }

FilterParameter& FilterParameter::Value(const string& value)
{ Attribute("Value", value); return *this; }

// ---------------------------------
// FilterParameters implementation
// ---------------------------------

FilterParameters::FilterParameters(void)
    : DataSetListElement<FilterParameter>("Parameters")
{ }

void FilterParameters::AddParameter(const FilterParameter& param)
{ AddChild(param); }

void FilterParameters::RemoveParameter(const FilterParameter& param)
{ RemoveChild(param); }

// --------------------------------
// Filter implementation
// --------------------------------

Filter::Filter(void)
    : DataSetElement("Filter")
{ AddChild(FilterParameters()); }

void Filter::AddParameter(const FilterParameter& param)
{ FilterParameterList().AddParameter(param); }

const FilterParameters& Filter::FilterParameterList(void) const
{
    try {
        return Child<FilterParameters>("Parameters");
    } catch (std::exception&) {
        return internal::NullParams();
    }
}

FilterParameters& Filter::FilterParameterList(void)
{
    if (!HasChild("Parameters"))
        AddChild(internal::NullParams());
    return Child<FilterParameters>("Parameters");
}

size_t Filter::NumFilterParameters(void) const
{ return FilterParameterList().NumChildren(); }

void Filter::RemoveParameter(const FilterParameter& param)
{ FilterParameterList().RemoveParameter(param); }

// --------------------------------
// Filters implementation
// --------------------------------

Filters::Filters(void)
    : DataSetListElement<Filter>("Filters")
{ }

void Filters::AddFilter(const Filter& filter)
{ AddChild(filter); }

size_t Filters::NumFilters(void) const
{ return NumChildren(); }

void Filters::RemoveFilter(const Filter& filter)
{ RemoveChild(filter); }
