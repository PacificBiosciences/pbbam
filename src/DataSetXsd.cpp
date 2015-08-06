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

#include "pbbam/DataSetXsd.h"
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

static map<XsdType, NamespaceInfo> DefaultRegistry(void)
{
    map<XsdType, NamespaceInfo> result;
    result[XsdType::NONE]                   = NamespaceInfo();
    result[XsdType::AUTOMATION_CONSTRAINTS] = NamespaceInfo("",       "http://pacificbiosciences.com/PacBioAutomationConstraints.xsd");
    result[XsdType::BASE_DATA_MODEL]        = NamespaceInfo("pbbase", "http://pacificbiosciences.com/PacBioBaseDataModel.xsd");
    result[XsdType::COLLECTION_METADATA]    = NamespaceInfo("pbmeta", "http://pacificbiosciences.com/PacBioCollectionMetadata.xsd");
    result[XsdType::COMMON_MESSAGES]        = NamespaceInfo("",       "http://pacificbiosciences.com/PacBioCommonMessages.xsd");
    result[XsdType::DATA_MODEL]             = NamespaceInfo("pbdm",   "http://pacificbiosciences.com/PacBioDataModel.xsd");
    result[XsdType::DATA_STORE]             = NamespaceInfo("",       "http://pacificbiosciences.com/PacBioDataStore.xsd");
    result[XsdType::DATASETS]               = NamespaceInfo("pbds",   "http://pacificbiosciences.com/PacBioDatasets.xsd");
    result[XsdType::DECL_DATA]              = NamespaceInfo("",       "http://pacificbiosciences.com/PacBioDeclData.xsd");
    result[XsdType::PART_NUMBERS]           = NamespaceInfo("pbpn",   "http://pacificbiosciences.com/PacBioPartNumbers.xsd");
    result[XsdType::PRIMARY_METRICS]        = NamespaceInfo("",       "http://pacificbiosciences.com/PacBioPrimaryMetrics.xsd");
    result[XsdType::REAGENT_KIT]            = NamespaceInfo("pbrk",   "http://pacificbiosciences.com/PacBioReagentKit.xsd");
    result[XsdType::RIGHTS_AND_ROLES]       = NamespaceInfo("",       "http://pacificbiosciences.com/PacBioRightsAndRoles.xsd");
    result[XsdType::SAMPLE_INFO]            = NamespaceInfo("pbsample", "http://pacificbiosciences.com/PacBioSampleInfo.xsd");
    result[XsdType::SEEDING_DATA]           = NamespaceInfo("",       "http://pacificbiosciences.com/PacBioSeedingData.xsd");
    return result;
}

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
    : data_(std::move(internal::DefaultRegistry()))
    , defaultXsdType_(XsdType::DATASETS)
{ }

NamespaceRegistry::NamespaceRegistry(const NamespaceRegistry &other)
    : data_(other.data_)
    , defaultXsdType_(other.defaultXsdType_)
{ }

NamespaceRegistry& NamespaceRegistry::operator=(const NamespaceRegistry& other)
{
    data_ = other.data_;
    defaultXsdType_ = other.defaultXsdType_;
    return *this;
}

NamespaceRegistry::~NamespaceRegistry(void) { }

const NamespaceInfo& NamespaceRegistry::DefaultNamespace(void) const
{ return Namespace(DefaultXsd()); }

XsdType NamespaceRegistry::DefaultXsd(void) const
{ return defaultXsdType_; }

const NamespaceInfo& NamespaceRegistry::Namespace(const XsdType& xsd) const
{ return data_.at(xsd); }

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

void NamespaceRegistry::Register(const XsdType& xsd, const NamespaceInfo& namespaceInfo)
{ data_[xsd] = namespaceInfo; }

void NamespaceRegistry::SetDefaultXsd(const XsdType& xsd)
{ defaultXsdType_ = xsd; }
