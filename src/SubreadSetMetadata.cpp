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

#include "pbbam/dataset/SubreadSetMetadata.h"
#include <exception>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

// empty, "null" components

static
const BioSampleReferencesMetadata& NullBioSampleReferences(void)
{
    static const BioSampleReferencesMetadata empty;
    return empty;
}

static
const BioSamplesMetadata& NullBioSamples(void)
{
    static const BioSamplesMetadata empty;
    return empty;
}

static
const CollectionsMetadata& NullCollections(void)
{
    static const CollectionsMetadata empty;
    return empty;
}

static
const CopyFilesMetadata& NullCopyFiles(void)
{
    static const CopyFilesMetadata empty;
    return empty;
}

static
const PrimaryMetadata& NullPrimary(void)
{
    static const PrimaryMetadata empty;
    return empty;
}

static
const RunDetailsMetadata& NullRunDetails(void)
{
    static const RunDetailsMetadata empty;
    return empty;
}

static
const WellSampleMetadata& NullWellSample(void)
{
    static const WellSampleMetadata empty;
    return empty;
}

} // namespace internal
} // namespace BAM
} // namespace PacBio

// -----------------------------------
// BioSampleReferenceMetadata implementation
// -----------------------------------

BioSampleReferencesMetadata::BioSampleReferencesMetadata(void)
    : internal::DataSetElement("BioSampleReferences")
{ }

// -----------------------------------
// BioSampleMetadata implementation
// -----------------------------------

BioSampleMetadata::BioSampleMetadata(void)
    : internal::DataSetElement("BioSample")
{ }

const string& BioSampleMetadata::CreatedAt(void) const
{ return FetchChildText("CreatedAt"); }

BioSampleMetadata& BioSampleMetadata::CreatedAt(const string& createdAt)
{ SetChildText("CreatedAt", createdAt); return *this; }

const string& BioSampleMetadata::UniqueId(void) const
{ return FetchChildText("UniqueId"); }

BioSampleMetadata& BioSampleMetadata::UniqueId(const string& uuid)
{ SetChildText("UniqueId", uuid); return *this; }

// -----------------------------------
// BioSamplesMetadata implementation
// -----------------------------------

BioSamplesMetadata::BioSamplesMetadata(void)
    : internal::DataSetListElement<BioSampleMetadata>("BioSamples")
{ }

BioSamplesMetadata& BioSamplesMetadata::AddBioSample(const BioSampleMetadata& bioSample)
{ AddChild(bioSample); return *this; }

BioSamplesMetadata& BioSamplesMetadata::RemoveBioSample(const BioSampleMetadata& bioSample)
{ RemoveChild(bioSample); return *this; }

// -----------------------------------
// CollectionMetadata implementation
// -----------------------------------

CollectionMetadata::CollectionMetadata(void)
    : internal::DataSetElement("Collection")
{ }

const string& CollectionMetadata::AutomationName(void) const
{ return FetchChildText("AutomationName"); }

CollectionMetadata& CollectionMetadata::AutomationName(const string& name)
{ SetChildText("AutomationName", name); return *this; }

const string& CollectionMetadata::CellIndex(void) const
{ return FetchChildText("CellIndex"); }

CollectionMetadata& CollectionMetadata::CellIndex(const string& index)
{ SetChildText("CellIndex", index); return *this; }

const string& CollectionMetadata::CellPac(void) const
{ return FetchChildText("CellPac"); }

CollectionMetadata& CollectionMetadata::CellPac(const string& pac)
{ SetChildText("CellPac", pac); return *this; }

const string& CollectionMetadata::Context(void) const
{ return FetchChildText("Context"); }

CollectionMetadata& CollectionMetadata::Context(const string& context)
{ SetChildText("Context", context); return *this; }

const string& CollectionMetadata::InstrCtrlVer(void) const
{ return FetchChildText("InstrCtrlVer"); }

CollectionMetadata& CollectionMetadata::InstrCtrlVer(const string& ver)
{ SetChildText("InstrCtrlVer", ver); return *this; }

const string& CollectionMetadata::InstrumentId(void) const
{ return FetchChildText("InstrumentId"); }

CollectionMetadata& CollectionMetadata::InstrumentId(const string& id)
{ SetChildText("InstrumentId", id); return *this; }

const string& CollectionMetadata::InstrumentName(void) const
{ return FetchChildText("InstrumentName"); }

CollectionMetadata& CollectionMetadata::InstrumentName(const string& name)
{ SetChildText("InstrumentName", name); return *this; }

const PrimaryMetadata& CollectionMetadata::Primary(void) const
{
    try {
        return Child<PrimaryMetadata>("Primary");
    } catch (std::exception&) {
        return internal::NullPrimary();
    }
}

PrimaryMetadata& CollectionMetadata::Primary(void)
{
    if (!HasChild("Primary"))
        AddChild(internal::NullPrimary());
    return Child<PrimaryMetadata>("Primary");
}

const RunDetailsMetadata& CollectionMetadata::RunDetails(void) const
{
    try {
        return Child<RunDetailsMetadata>("RunDetails");
    } catch (std::exception&) {
        return internal::NullRunDetails();
    }
}

RunDetailsMetadata& CollectionMetadata::RunDetails(void)
{
    if (!HasChild("RunDetails"))
        AddChild(internal::NullRunDetails());
    return Child<RunDetailsMetadata>("RunDetails");
}

const string& CollectionMetadata::SigProcVer(void) const
{ return FetchChildText("SigProcVer"); }

CollectionMetadata& CollectionMetadata::SigProcVer(const string& ver)
{ SetChildText("SigProcVer", ver); return *this; }

const WellSampleMetadata& CollectionMetadata::WellSample(void) const
{
    try {
        return Child<WellSampleMetadata>("WellSample");
    } catch (std::exception&) {
        return internal::NullWellSample();
    }
}

WellSampleMetadata& CollectionMetadata::WellSample(void)
{
    if (!HasChild("WellSample"))
        AddChild(internal::NullWellSample());
    return Child<WellSampleMetadata>("WellSample");
}

// -----------------------------------
// CollectionsMetadata implementation
// -----------------------------------

CollectionsMetadata::CollectionsMetadata(void)
    : internal::DataSetListElement<CollectionMetadata>("Collections")
{ }

CollectionsMetadata& CollectionsMetadata::AddCollection(const CollectionMetadata& collection)
{ AddChild(collection); return *this; }

CollectionsMetadata& CollectionsMetadata::RemoveCollection(const CollectionMetadata& collection)
{ RemoveChild(collection); return *this; }

// -----------------------------------
// CopyFilesMetadata implementation
// -----------------------------------

 CopyFilesMetadata::CopyFilesMetadata(void)
     : internal::DataSetElement("CopyFiles")
 { }

// -----------------------------------
// PrimaryMetadata implementation
// -----------------------------------

PrimaryMetadata::PrimaryMetadata(void)
    : internal::DataSetElement("Primary")
{ }

const string& PrimaryMetadata::AutomationName(void) const
{ return FetchChildText("AutomationName"); }

PrimaryMetadata& PrimaryMetadata::AutomationName(const string& name)
{ SetChildText("AutomationName", name); return *this; }

const string& PrimaryMetadata::CollectionPathUri(void) const
{ return FetchChildText("CollectionPathUri"); }

PrimaryMetadata& PrimaryMetadata::CollectionPathUri(const string& uri)
{ SetChildText("CollectionPathUri", uri); return *this; }

const string& PrimaryMetadata::ContigFileName(void) const
{ return FetchChildText("ContigFileName"); }

PrimaryMetadata& PrimaryMetadata::ContigFileName(const string& name)
{ SetChildText("ContigFileName", name); return *this; }

const CopyFilesMetadata& PrimaryMetadata::CopyFiles(void) const
{
    try {
        return Child<CopyFilesMetadata>("CopyFiles");
    } catch (std::exception&) {
        return internal::NullCopyFiles();
    }
}

CopyFilesMetadata& PrimaryMetadata::CopyFiles(void)
{
    if (!HasChild("CopyFiles"))
        AddChild(internal::NullCopyFiles());
    return Child<CopyFilesMetadata>("CopyFiles");
}

const string& PrimaryMetadata::ResultsFolder(void) const
{ return FetchChildText("ResultsFolder"); }

PrimaryMetadata& PrimaryMetadata::ResultsFolder(const string& folder)
{ SetChildText("ResultsFolder", folder); return *this; }

const string& PrimaryMetadata::SequencingCondition(void) const
{ return FetchChildText("SequencingCondition"); }

PrimaryMetadata& PrimaryMetadata::SequencingCondition(const string& condition)
{ SetChildText("SequencingCondition", condition); return *this; }

// -----------------------------------
// RunDetailsMetadata implementation
// -----------------------------------

RunDetailsMetadata::RunDetailsMetadata(void)
    : internal::DataSetElement("RunDetails")
{ }

const string& RunDetailsMetadata::Name(void) const
{ return FetchChildText("Name"); }

RunDetailsMetadata& RunDetailsMetadata::Name(const string& name)
{ SetChildText("Name", name); return *this; }

const string& RunDetailsMetadata::RunId(void) const
{ return FetchChildText("RunId"); }

RunDetailsMetadata& RunDetailsMetadata::RunId(const string& runId)
{ SetChildText("RunId", runId); return *this; }

// -----------------------------------
// SubreadSetMetadata implementation
// -----------------------------------

SubreadSetMetadata::SubreadSetMetadata(void)
    : DataSetMetadataBase()
{ }

const BioSamplesMetadata& SubreadSetMetadata::BioSamples(void) const
{
    try {
        return Child<BioSamplesMetadata>("BioSamples");
    } catch (std::exception&) {
        return internal::NullBioSamples();
    }
}

BioSamplesMetadata& SubreadSetMetadata::BioSamples(void)
{
    if (!HasChild("BioSamples"))
        AddChild(internal::NullBioSamples());
    return Child<BioSamplesMetadata>("BioSamples");
}

const CollectionsMetadata& SubreadSetMetadata::Collections(void) const
{
    try {
        return Child<CollectionsMetadata>("Collections");
    } catch (std::exception&) {
        return internal::NullCollections();
    }
}

CollectionsMetadata& SubreadSetMetadata::Collections(void)
{
    if (!HasChild("Collections"))
        AddChild(internal::NullCollections());
    return Child<CollectionsMetadata>("Collections");
}

// -----------------------------------
// WellSampleMetadata implementation
// -----------------------------------

WellSampleMetadata::WellSampleMetadata(void)
    : internal::DataSetElement("WellSample")
{ }

const BioSampleReferencesMetadata& WellSampleMetadata::BioSampleReferences(void) const
{
    try {
        return Child<BioSampleReferencesMetadata>("BioSampleReferences");
    } catch (std::exception&) {
        return internal::NullBioSampleReferences();
    }
}

BioSampleReferencesMetadata& WellSampleMetadata::BioSampleReferences(void)
{
    if (!HasChild("BioSampleReferences"))
        AddChild(internal::NullBioSampleReferences());
    return Child<BioSampleReferencesMetadata>("BioSampleReferences");
}

const string& WellSampleMetadata::Comments(void) const
{ return FetchChildText("Comments"); }

WellSampleMetadata& WellSampleMetadata::Comments(const string& comments)
{ SetChildText("Comments", comments); return *this; }

const string& WellSampleMetadata::Concentration(void) const
{ return FetchChildText("Concentration"); }

WellSampleMetadata& WellSampleMetadata::Concentration(const string& concentration)
{ SetChildText("Concentration", concentration); return *this; }

const string& WellSampleMetadata::PlateId(void) const
{ return FetchChildText("PlateId"); }

WellSampleMetadata& WellSampleMetadata::PlateId(const string& plateId)
{ SetChildText("PlateId", plateId); return *this; }

const string& WellSampleMetadata::SampleReuseEnabled(void) const
{ return FetchChildText("SampleReuseEnabled"); }

WellSampleMetadata& WellSampleMetadata::SampleReuseEnabled(const string& enabled)
{ SetChildText("SampleReuseEnabled", enabled); return *this; }

const string& WellSampleMetadata::SizeSelectionEnabled(void) const
{ return FetchChildText("SizeSelectionEnabled"); }

WellSampleMetadata& WellSampleMetadata::SizeSelectionEnabled(const string& enabled)
{ SetChildText("SizeSelectionEnabled", enabled); return *this; }

const string& WellSampleMetadata::StageHotstartEnabled(void) const
{ return FetchChildText("StageHotstartEnabled"); }

WellSampleMetadata& WellSampleMetadata::StageHotstartEnabled(const string& enabled)
{ SetChildText("StageHotstartEnabled", enabled); return *this; }

const string& WellSampleMetadata::UniqueId(void) const
{ return FetchChildText("UniqueId"); }

WellSampleMetadata& WellSampleMetadata::UniqueId(const string& uuid)
{ SetChildText("UniqueId", uuid); return *this; }

const string& WellSampleMetadata::UseCount(void) const
{ return FetchChildText("UseCount"); }

WellSampleMetadata& WellSampleMetadata::UseCount(const string& count)
{ SetChildText("UseCount", count); return *this; }

const string& WellSampleMetadata::WellName(void) const
{ return FetchChildText("WellName"); }

WellSampleMetadata& WellSampleMetadata::WellName(const string& name)
{ SetChildText("WellName", name); return *this; }
