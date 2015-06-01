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

#include "pbbam/dataset/SubDataSets.h"
using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

// empty, "null" lists
static const Filters NullFilters = Filters();

} // namespace internal
} // namespace BAM
} // namespace PacBio

// ---------------------------
// SubDataSet implementation
// ---------------------------

SubDataSet::SubDataSet(void)
    : DataSetElement("DataSet")
{  AddChild(Filters("Filters")); }

void SubDataSet::AddFilter(const Filter& filter)
{ FilterList().AddFilter(filter); }

string SubDataSet::CreatedAt(void) const
{ return Attribute("CreatedAt"); }

SubDataSet& SubDataSet::CreatedAt(const string& timestamp)
{ Attribute("CreatedAt", timestamp); return *this; }

const Filters& SubDataSet::FilterList(void) const
{
    try {
        return Child<Filters>("Filters");
    } catch (std::exception&) {
        return internal::NullFilters;
    }
}

Filters& SubDataSet::FilterList(void)
{
    if (!HasChild("Filters"))
        AddChild(internal::NullFilters);
    return Child<Filters>("Filters");
}

std::string SubDataSet::Name(void) const
{ return Attribute("Name"); }

SubDataSet& SubDataSet::Name(const std::string& name)
{ Attribute("Name", name); return *this; }

void SubDataSet::RemoveFilter(const Filter& filter)
{ FilterList().RemoveFilter(filter); }

std::string SubDataSet::Tags(void) const
{ return Attribute("Tags"); }

SubDataSet& SubDataSet::Tags(const std::string& tags)
{ Attribute("Tags", tags); return *this; }

std::string SubDataSet::UniqueId(void) const
{ return Attribute("UniqueId"); }

SubDataSet& SubDataSet::UniqueId(const std::string& uuid)
{ Attribute("UniqueId", uuid); return *this;}

std::string SubDataSet::Version(void) const
{ return Attribute("Version"); }

SubDataSet& SubDataSet::Version(const std::string& version)
{ Attribute("Version", version); return *this;}

// ----------------------------
// SubDataSets implementation
// ----------------------------

SubDataSets::SubDataSets(void)
    : DataSetListElement<SubDataSet>("DataSets")
{ }

void SubDataSets::AddSubDataSet(const SubDataSet& subdataset)
{ AddChild(subdataset); }

void SubDataSets::RemoveSubDataSet(const SubDataSet& subdataset)
{ RemoveChild(subdataset); }
