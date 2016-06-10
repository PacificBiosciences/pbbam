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
/// \file DataSetXsd.cpp
/// \brief Implements the XSD- and namespace-related classes for DataSetXML.
//
// Author: Derek Barnett

#include "pbbam/DataSetXsd.h"
#include <unordered_map>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

static map<XsdType, NamespaceInfo> DefaultRegistry(void)
{
    const auto result = map<XsdType, NamespaceInfo>
    {
        { XsdType::NONE,                   NamespaceInfo{ "", "" } },
        { XsdType::AUTOMATION_CONSTRAINTS, NamespaceInfo{ "",       "http://pacificbiosciences.com/PacBioAutomationConstraints.xsd" } },
        { XsdType::BASE_DATA_MODEL,        NamespaceInfo{ "pbbase", "http://pacificbiosciences.com/PacBioBaseDataModel.xsd" } },
        { XsdType::COLLECTION_METADATA,    NamespaceInfo{ "pbmeta", "http://pacificbiosciences.com/PacBioCollectionMetadata.xsd" } },
        { XsdType::COMMON_MESSAGES,        NamespaceInfo{ "",       "http://pacificbiosciences.com/PacBioCommonMessages.xsd" } },
        { XsdType::DATA_MODEL,             NamespaceInfo{ "pbdm",   "http://pacificbiosciences.com/PacBioDataModel.xsd" } },
        { XsdType::DATA_STORE,             NamespaceInfo{ "",       "http://pacificbiosciences.com/PacBioDataStore.xsd" } },
        { XsdType::DATASETS,               NamespaceInfo{ "pbds",   "http://pacificbiosciences.com/PacBioDatasets.xsd" } },
        { XsdType::DECL_DATA,              NamespaceInfo{ "",       "http://pacificbiosciences.com/PacBioDeclData.xsd" } },
        { XsdType::PART_NUMBERS,           NamespaceInfo{ "pbpn",   "http://pacificbiosciences.com/PacBioPartNumbers.xsd" } },
        { XsdType::PRIMARY_METRICS,        NamespaceInfo{ "",       "http://pacificbiosciences.com/PacBioPrimaryMetrics.xsd" } },
        { XsdType::REAGENT_KIT,            NamespaceInfo{ "pbrk",   "http://pacificbiosciences.com/PacBioReagentKit.xsd" } },
        { XsdType::RIGHTS_AND_ROLES,       NamespaceInfo{ "",       "http://pacificbiosciences.com/PacBioRightsAndRoles.xsd" } },
        { XsdType::SAMPLE_INFO,            NamespaceInfo{ "pbsample", "http://pacificbiosciences.com/PacBioSampleInfo.xsd" } },
        { XsdType::SEEDING_DATA,           NamespaceInfo{ "",       "http://pacificbiosciences.com/PacBioSeedingData.xsd" } }
    };
    return result;
}

static const auto elementRegistry = unordered_map<string, XsdType>
{
    // 'pbbase' elements
    //
    { "AutomationParameter" ,  XsdType::BASE_DATA_MODEL },
    { "AutomationParameters" , XsdType::BASE_DATA_MODEL },
    { "BinCount" ,             XsdType::BASE_DATA_MODEL },
    { "BinCounts" ,            XsdType::BASE_DATA_MODEL },
    { "BinLabel" ,             XsdType::BASE_DATA_MODEL },
    { "BinLabels" ,            XsdType::BASE_DATA_MODEL },
    { "BinWidth" ,             XsdType::BASE_DATA_MODEL },
    { "ExternalResource" ,     XsdType::BASE_DATA_MODEL },
    { "ExternalResources" ,    XsdType::BASE_DATA_MODEL },
    { "FileIndex" ,            XsdType::BASE_DATA_MODEL },
    { "FileIndices" ,          XsdType::BASE_DATA_MODEL },
    { "MaxBinValue" ,          XsdType::BASE_DATA_MODEL },
    { "MaxOutlierValue" ,      XsdType::BASE_DATA_MODEL },
    { "MetricDescription" ,    XsdType::BASE_DATA_MODEL },
    { "NumBins" ,              XsdType::BASE_DATA_MODEL },
    { "Properties" ,           XsdType::BASE_DATA_MODEL },
    { "Property" ,             XsdType::BASE_DATA_MODEL },
    { "Sample95thPct" ,        XsdType::BASE_DATA_MODEL },
    { "SampleMean" ,           XsdType::BASE_DATA_MODEL },
    { "SampleMed" ,            XsdType::BASE_DATA_MODEL },
    { "SampleSize" ,           XsdType::BASE_DATA_MODEL },
    { "SampleStd" ,            XsdType::BASE_DATA_MODEL },

    // 'pbds' elements
    //
    { "AdapterDimerFraction",  XsdType::DATASETS },
    { "AlignmentSet",          XsdType::DATASETS },
    { "BarcodeConstruction",   XsdType::DATASETS },
    { "BarcodeSet",            XsdType::DATASETS },
    { "ConsensusAlignmentSet", XsdType::DATASETS },
    { "ConsensusReadSet",      XsdType::DATASETS },
    { "Contig",                XsdType::DATASETS },
    { "Contigs",               XsdType::DATASETS },
    { "ContigSet",             XsdType::DATASETS },
    { "ControlReadLenDist",    XsdType::DATASETS },
    { "ControlReadQualDist",   XsdType::DATASETS },
    { "DataSetMetdata",        XsdType::DATASETS },
    { "DataSet",               XsdType::DATASETS },
    { "DataSets",              XsdType::DATASETS },
    { "Filter",                XsdType::DATASETS },
    { "Filters",               XsdType::DATASETS },
    { "HdfSubreadSet",         XsdType::DATASETS },
    { "InsertReadLenDist",     XsdType::DATASETS },
    { "InsertReadQualDist" ,   XsdType::DATASETS },
    { "MedianInsertDist",      XsdType::DATASETS },
    { "NumRecords",            XsdType::DATASETS },
    { "NumSequencingZmws",     XsdType::DATASETS },
    { "Organism",              XsdType::DATASETS },
    { "ParentTool",            XsdType::DATASETS },
    { "Ploidy",                XsdType::DATASETS },
    { "ProdDist",              XsdType::DATASETS },
    { "Provenance",            XsdType::DATASETS },
    { "ReadLenDist",           XsdType::DATASETS },
    { "ReadQualDist",          XsdType::DATASETS },
    { "ReadTypeDist",          XsdType::DATASETS },
    { "ReferenceSet",          XsdType::DATASETS },
    { "ShortInsertFraction",   XsdType::DATASETS },
    { "SubreadSet",            XsdType::DATASETS },
    { "SummaryStats",          XsdType::DATASETS },
    { "TotalLength",           XsdType::DATASETS },

    // 'pbmeta' elements
    //
    { "Automation",           XsdType::COLLECTION_METADATA },
    { "AutomationName",       XsdType::COLLECTION_METADATA },
    { "CellIndex",            XsdType::COLLECTION_METADATA },
    { "CellPac",              XsdType::COLLECTION_METADATA },
    { "CollectionFileCopy",   XsdType::COLLECTION_METADATA },
    { "CollectionMetadata",   XsdType::COLLECTION_METADATA },
    { "CollectionNumber",     XsdType::COLLECTION_METADATA },
    { "CollectionPathUri",    XsdType::COLLECTION_METADATA },
    { "Collections",          XsdType::COLLECTION_METADATA },
    { "Concentration",        XsdType::COLLECTION_METADATA },
    { "ConfigFileName",       XsdType::COLLECTION_METADATA },
    { "CopyFiles",            XsdType::COLLECTION_METADATA },
    { "InstCtrlVer",          XsdType::COLLECTION_METADATA },
    { "MetricsVerbosity",     XsdType::COLLECTION_METADATA },
    { "Name",                 XsdType::COLLECTION_METADATA },
    { "OutputOptions",        XsdType::COLLECTION_METADATA },
    { "PlateId",              XsdType::COLLECTION_METADATA },
    { "Primary",              XsdType::COLLECTION_METADATA },
    { "Readout",              XsdType::COLLECTION_METADATA },
    { "ResultsFolder",        XsdType::COLLECTION_METADATA },
    { "RunDetails",           XsdType::COLLECTION_METADATA },
    { "RunId",                XsdType::COLLECTION_METADATA },
    { "SampleReuseEnabled",   XsdType::COLLECTION_METADATA },
    { "SequencingCondition",  XsdType::COLLECTION_METADATA },
    { "SigProcVer",           XsdType::COLLECTION_METADATA },
    { "SizeSelectionEnabled", XsdType::COLLECTION_METADATA },
    { "StageHotstartEnabled", XsdType::COLLECTION_METADATA },
    { "UseCount",             XsdType::COLLECTION_METADATA },
    { "WellName",             XsdType::COLLECTION_METADATA },
    { "WellSample",           XsdType::COLLECTION_METADATA },

    // 'pbsample' elements
    //
    { "BioSample",         XsdType::SAMPLE_INFO },
    { "BioSamplePointer",  XsdType::SAMPLE_INFO },
    { "BioSamplePointers", XsdType::SAMPLE_INFO },
    { "BioSamples",        XsdType::SAMPLE_INFO }
};

} // namespace internal
} // namespace BAM
} // namespace PacBio

// ---------------
// NamespaceInfo
// ---------------

NamespaceInfo::NamespaceInfo(void) { }

NamespaceInfo::NamespaceInfo(const string& name,
                             const string& uri)
    : name_(name)
    , uri_(uri)
{ }

// -------------------
// NamespaceRegistry
// -------------------

NamespaceRegistry::NamespaceRegistry(void)
    : data_(internal::DefaultRegistry())
    , defaultXsdType_(XsdType::DATASETS)
{ }

NamespaceRegistry::NamespaceRegistry(const NamespaceRegistry &other)
    : data_(other.data_)
    , defaultXsdType_(other.defaultXsdType_)
{ }

NamespaceRegistry::NamespaceRegistry(NamespaceRegistry &&other)
    : data_(std::move(other.data_))
    , defaultXsdType_(std::move(other.defaultXsdType_))
{ }

NamespaceRegistry& NamespaceRegistry::operator=(const NamespaceRegistry& other)
{
    data_ = other.data_;
    defaultXsdType_ = other.defaultXsdType_;
    return *this;
}

NamespaceRegistry& NamespaceRegistry::operator=(NamespaceRegistry&& other)
{
    data_ = std::move(other.data_);
    defaultXsdType_ = std::move(other.defaultXsdType_);
    return *this;
}

NamespaceRegistry::~NamespaceRegistry(void) { }

const NamespaceInfo& NamespaceRegistry::DefaultNamespace(void) const
{ return Namespace(DefaultXsd()); }

XsdType NamespaceRegistry::DefaultXsd(void) const
{ return defaultXsdType_; }

const NamespaceInfo& NamespaceRegistry::Namespace(const XsdType& xsd) const
{ return data_.at(xsd); }

void NamespaceRegistry::Register(const XsdType& xsd, const NamespaceInfo& namespaceInfo)
{ data_[xsd] = namespaceInfo; }

void NamespaceRegistry::SetDefaultXsd(const XsdType& xsd)
{ defaultXsdType_ = xsd; }

XsdType NamespaceRegistry::XsdForElement(const std::string& elementLabel) const
{
    const auto iter = internal::elementRegistry.find(elementLabel);
    return (iter == internal::elementRegistry.cend() ? XsdType::NONE : iter->second);
}

XsdType NamespaceRegistry::XsdForUri(const std::string& uri) const
{
    map<XsdType, NamespaceInfo>::const_iterator iter = data_.cbegin();
    map<XsdType, NamespaceInfo>::const_iterator end  = data_.cend();
    for ( ; iter != end; ++iter ) {
        const NamespaceInfo& info = iter->second;
        if (info.Uri() == uri)
            return iter->first;
    }
    return XsdType::NONE;
}
