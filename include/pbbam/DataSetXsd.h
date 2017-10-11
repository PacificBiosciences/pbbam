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
/// \file DataSetXsd.h
/// \brief Defines the XSD- and namespace-related classes for DataSetXML.
//
// Author: Derek Barnett

#ifndef DATASETXSD_H
#define DATASETXSD_H

#include "pbbam/Config.h"
#include <map>
#include <string>

namespace PacBio {
namespace BAM {

/// \brief The XsdType enum defines the supported XSD namespaces.
///
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

/// \brief The NamespaceInfo class provides XML namespace info (prefix & URI).
///
class PBBAM_EXPORT NamespaceInfo
{
public:
    /// \brief Creates an empty entry.
    ///
    /// This constructor only exists for STL container compatibility.
    ///
    NamespaceInfo() = default;

    /// \brief Creates a valid info entry.
    NamespaceInfo(std::string name,
                  std::string uri);

public:
    /// \brief Fetches namespace name (i.e. prefix)
    const std::string& Name() const { return name_; }

    /// \brief Fetches namespace URI.
    const std::string& Uri() const { return uri_; }

private:
    std::string name_;
    std::string uri_;
};

/// \brief The NamespaceRegistry class provides a per-dataset registry of XML
///        namespace information.
///
/// This is used to format XML output - properly prefixing element labels with
/// namespace as appropriate.
///
class PBBAM_EXPORT NamespaceRegistry
{
public:
    /// \name Constructors & Related Methods
    /// \{

    NamespaceRegistry();
    NamespaceRegistry(const NamespaceRegistry& other) = default;
    NamespaceRegistry(NamespaceRegistry&& other) = default;
    NamespaceRegistry& operator=(const NamespaceRegistry& other) = default;
    NamespaceRegistry& operator=(NamespaceRegistry&& other) = default;
    ~NamespaceRegistry() = default;

    /// \}

public:
    /// \name Registry Access
    /// \{

    /// \brief Fetches namespace info for the dataset's default XSD type.
    const NamespaceInfo& DefaultNamespace() const;

    /// \brief Fetches dataset's default XSD type.
    XsdType DefaultXsd() const;

    /// \brief Fetches namespace info for the requested XSD type.
    const NamespaceInfo& Namespace(const XsdType& xsd) const;

    /// \brief Registers namespace info for a particular XSD type.
    void Register(const XsdType& xsd, const NamespaceInfo& namespaceInfo);

    /// \brief Updates dataset's default XSD type.
    void SetDefaultXsd(const XsdType& xsd);

    /// \brief Fetches the XSD type for \p elementLabel.
    XsdType XsdForElement(const std::string& elementLabel) const;

    /// \brief Fetches the XSD type for a particular URI.
    XsdType XsdForUri(const std::string& uri) const;

    /// \}

private:
    std::map<XsdType, NamespaceInfo> data_;
    XsdType defaultXsdType_ = XsdType::DATASETS;
};

} // namespace PacBio
} // namespace BAM

#endif // DATASETXSD_H
