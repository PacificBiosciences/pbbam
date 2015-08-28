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

#ifndef DATASETXSD_H
#define DATASETXSD_H

#include "pbbam/Config.h"
#include <map>
#include <string>

namespace PacBio {
namespace BAM {

enum class XsdType
{
    NONE

  , AUTOMATION_CONSTRAINTS
  , BASE_DATA_MODEL
  , COLLECTION_METADATA
  , COMMON_MESSAGES
  , DATA_MODEL
  , DATA_STORE
  , DATASETS
  , DECL_DATA
  , PART_NUMBERS
  , PRIMARY_METRICS
  , REAGENT_KIT
  , RIGHTS_AND_ROLES
  , SAMPLE_INFO
  , SEEDING_DATA
};

class PBBAM_EXPORT NamespaceInfo
{
public:
    NamespaceInfo(void);
    NamespaceInfo(const std::string& name,
                  const std::string& uri);

public:
    const std::string& Name(void) const { return name_; }
    const std::string& Uri(void) const { return uri_; }

private:
    std::string name_;
    std::string uri_;
};

class PBBAM_EXPORT NamespaceRegistry
{
public:
    NamespaceRegistry(void);
    NamespaceRegistry(const NamespaceRegistry& other);
    NamespaceRegistry(NamespaceRegistry&& other);
    NamespaceRegistry& operator=(const NamespaceRegistry& other);
    NamespaceRegistry& operator=(NamespaceRegistry&& other);
    ~NamespaceRegistry(void);

public:
    const NamespaceInfo& DefaultNamespace(void) const;
    XsdType DefaultXsd(void) const;
    const NamespaceInfo& Namespace(const XsdType& xsd) const;

    XsdType XsdForUri(const std::string& uri) const;

public:
    void Register(const XsdType& xsd, const NamespaceInfo& namespaceInfo);
    void SetDefaultXsd(const XsdType& xsd);

private:
    std::map<XsdType, NamespaceInfo> data_;
    XsdType defaultXsdType_;
};

} // namespace PacBio
} // namespace BAM

#endif // DATASETXSD_H
