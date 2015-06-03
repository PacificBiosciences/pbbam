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

#include "pbbam/dataset/DataSetBase.h"
#include "pbbam/dataset/DataSetMetadata.h"
#include "DataSetIO.h"
#include "StringUtils.h"
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <set>
#include <vector>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

// empty, "null" components
static
const Filters& NullFilters(void)
{
    static const Filters empty;
    return empty;
}

static
const ExternalDataReferences& NullRefs(void)
{
    static const ExternalDataReferences empty;
    return empty;
}

static
const SubDataSets& NullSubDataSets(void)
{
    static const SubDataSets empty;
    return empty;
}

} // namespace internal
} // namespace BAM
} // namespace PacBio

DataSetBase::DataSetBase(void)
    : DataSetElement("DataSet")
{ }

DataSetBase::DataSetBase(const DataSetType& type)
    : DataSetElement(NameForType(type))
{ }

DataSetBase::DataSetBase(const string& filename)
    : DataSetElement("DataSet")
{ *this = std::move(DataSetIO::FromUri(filename)); }

DataSetBase::DataSetBase(const vector<string>& uris)
    : DataSetElement("DataSet")
{ *this = std::move(DataSetIO::FromUris(uris)); }

DataSetBase::DataSetBase(const DataSetBase& other)
    : DataSetElement(other)
{ }

DataSetBase::DataSetBase(DataSetBase&& other)
    : DataSetElement(std::move(other))
{ }

DataSetBase& DataSetBase::operator=(const DataSetBase& other)
{ DataSetElement::operator=(other); return *this; }

DataSetBase& DataSetBase::operator=(DataSetBase&& other)
{ DataSetElement::operator=(std::move(other)); return *this; }

DataSetBase::~DataSetBase(void) { }

DataSetBase DataSetBase::operator+(const DataSetBase& other) const
{
    DataSetBase result;
    result += other;
    return result;
}

DataSetBase& DataSetBase::operator+=(const DataSetBase& other)
{
    // fail on conflicting filters, just for now to simplify
    const Filters& filters = FilterList();
    const Filters& otherFilters = other.FilterList();
    if (filters != otherFilters )
        throw std::exception();

    // only keep unique resource ids
    const ExternalDataReferences& myRefs = ExternalDataReferenceList();
    const ExternalDataReferences& otherRefs = other.ExternalDataReferenceList();

    set<std::string> myResourceIds;
    for (auto& ref : myRefs)
        myResourceIds.insert(ref.ResourceId());

    vector<size_t> newResourceIndices;
    const size_t numOtherResourceIds = otherRefs.Size();
    for (size_t i = 0; i < numOtherResourceIds; ++i) {
        const string& resourceId = otherRefs[i].ResourceId();
        auto found = myResourceIds.find(resourceId);
        if (found == myResourceIds.cend())
            newResourceIndices.push_back(i);
    }

    for (size_t index : newResourceIndices)
        AddExternalDataReference( otherRefs[index] );

    return *this;
}

DataSetBase& DataSetBase::AddExternalDataReference(const ExternalDataReference& ref)
{ ExternalDataReferenceList().AddExternalRef(ref); return *this; }

DataSetBase& DataSetBase::AddFilter(const Filter& filter)
{ FilterList().AddFilter(filter); return *this; }

DataSetBase &DataSetBase::AddSubDataSet(const SubDataSet& subdataset)
{ SubDataSetList().AddSubDataSet(subdataset); return *this; }

const string& DataSetBase::CreatedAt(void) const
{ return Attribute("CreatedAt"); }

DataSetBase& DataSetBase::CreatedAt(const string& timestamp)
{ Attribute("CreatedAt", timestamp); return *this; }

const ExternalDataReferences& DataSetBase::ExternalDataReferenceList(void) const
{
    try {
        return Child<ExternalDataReferences>("ExternalDataReferences");
    } catch (std::exception&) {
        return internal::NullRefs();
    }
}

ExternalDataReferences& DataSetBase::ExternalDataReferenceList(void)
{
    if (!HasChild("ExternalDataReferences"))
        AddChild(internal::NullRefs());
    return Child<ExternalDataReferences>("ExternalDataReferences");
}

const Filters& DataSetBase::FilterList(void) const
{
    try {
        return Child<Filters>("Filters");
    } catch (std::exception&) {
        return internal::NullFilters();
    }
}

Filters& DataSetBase::FilterList(void)
{
    if (!HasChild("Filters"))
        AddChild(internal::NullFilters());
    return Child<Filters>("Filters");
}

const string& DataSetBase::MetaType(void) const
{ return Attribute("MetaType"); }

DataSetBase& DataSetBase::MetaType(const string& metatype)
{ Attribute("MetaType", metatype); return *this; }

const string& DataSetBase::Name(void) const
{ return Attribute("Name"); }

DataSetBase& DataSetBase::Name(const string& name)
{ Attribute("Name", name); return *this; }

string DataSetBase::NameForType(const DataSetType type)
{
    switch (type) {
        case DataSetType::GENERIC      : return "DataSet";
        case DataSetType::ALIGNMENTSET : return "AlignmentSet";
        case DataSetType::BARCODESET   : return "BarcodeSet";
        case DataSetType::CCSREADSET   : return "CCSreadSet";  // TODO: is this correct caps?
        case DataSetType::CONTIGSET    : return "ContigSet";
        case DataSetType::REFERENCESET : return "ReferenceSet";
        case DataSetType::SUBREADSET   : return "SubreadSet";
        default:
            throw std::exception();
    }
    // unreachable
    return string();
}

size_t DataSetBase::NumExternalDataReferences(void) const
{
    if (!HasChild("ExternalDataReferences"))
        return 0;
    return ExternalDataReferenceList().Size();
}

size_t DataSetBase::NumFilters(void) const
{
    if (!HasChild("Filters"))
        return 0;
    return FilterList().Size();
}

size_t DataSetBase::NumSubDataSets(void) const
{
    if (!HasChild("DataSets"))
        return 0;
    return SubDataSetList().Size();
}

DataSetBase& DataSetBase::RemoveExternalDataReference(const ExternalDataReference& ref)
{ ExternalDataReferenceList().RemoveExternalRef(ref); return *this; }

DataSetBase& DataSetBase::RemoveFilter(const Filter& filter)
{ FilterList().RemoveFilter(filter); return *this; }

DataSetBase& DataSetBase::RemoveSubDataSet(const SubDataSet& subdataset)
{ SubDataSetList().RemoveSubDataSet(subdataset); return *this; }

const SubDataSets& DataSetBase::SubDataSetList(void) const
{
    try {
        return Child<SubDataSets>("DataSets");
    } catch (std::exception&) {
        return internal::NullSubDataSets();
    }
}

SubDataSets& DataSetBase::SubDataSetList(void)
{
    if (!HasChild("DataSets"))
        AddChild(internal::NullSubDataSets());
    return Child<SubDataSets>("DataSets");
}

const string& DataSetBase::Tags(void) const
{ return Attribute("Tags"); }

DataSetBase& DataSetBase::Tags(const string& tags)
{ Attribute("Tags", tags); return *this; }

const DataSetType DataSetBase::Type(void) const
{ return TypeForName(Label()); }

DataSetBase& DataSetBase::Type(DataSetType type)
{ Label(NameForType(type)); return *this; }

DataSetType DataSetBase::TypeForName(const string& name)
{
    static map<string, DataSetType> lookup;
    if (lookup.empty()) {
        lookup["AlignmentSet"] = DataSetType::ALIGNMENTSET;
        lookup["BarcodeSet"]   = DataSetType::BARCODESET;
        lookup["CCSreadSet"]   = DataSetType::CCSREADSET; // TODO: is this correct caps?
        lookup["ContigSet"]    = DataSetType::CONTIGSET;
        lookup["DataSet"]      = DataSetType::GENERIC;
        lookup["ReferenceSet"] = DataSetType::REFERENCESET;
        lookup["SubreadSet"]   = DataSetType::SUBREADSET;
    }
    auto found = lookup.find(name);
    if (found == lookup.cend())
        throw std::exception();
    return found->second;
}

const string& DataSetBase::UniqueId(void) const
{ return Attribute("UniqueId"); }

DataSetBase& DataSetBase::UniqueId(const string& uuid)
{ Attribute("UniqueId", uuid); ; return *this; }

const string& DataSetBase::Version(void) const
{ return Attribute("Version"); }

DataSetBase& DataSetBase::Version(const string& version)
{ Attribute("Version", version); return *this; }

void DataSetBase::Write(const std::string& fn) const
{ DataSetIO::ToFile(*this, fn); }

void DataSetBase::WriteToStderr(void) const
{ WriteToStream(std::cerr); }

void DataSetBase::WriteToStdout(void) const
{ WriteToStream(std::cout); }

void DataSetBase::WriteToStream(std::ostream& out) const
{ DataSetIO::ToStream(*this, out); }
