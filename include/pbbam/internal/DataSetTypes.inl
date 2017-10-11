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
/// \file DataSetTypes.inl
/// \brief Inline implementations for the public DataSet component classes.
//
// Author: Derek Barnett

#include "pbbam/DataSetTypes.h"

namespace PacBio {
namespace BAM {

// -------------
// DataSetBase
// --------------

inline const NamespaceRegistry& DataSetBase::Namespaces() const
{ return registry_; }

inline NamespaceRegistry& DataSetBase::Namespaces()
{ return registry_; }

// ---------------------
// DataSetMetadata
// ---------------------

inline const std::string& DataSetMetadata::NumRecords() const
{ return ChildText("NumRecords"); }

inline std::string& DataSetMetadata::NumRecords()
{ return ChildText("NumRecords"); }

inline DataSetMetadata& DataSetMetadata::NumRecords(const std::string& numRecords)
{ ChildText("NumRecords", numRecords); return *this; }

inline const std::string& DataSetMetadata::TotalLength() const
{ return ChildText("TotalLength"); }

inline std::string& DataSetMetadata::TotalLength()
{ return ChildText("TotalLength"); }

inline DataSetMetadata& DataSetMetadata::TotalLength(const std::string& totalLength)
{ ChildText("TotalLength", totalLength); return *this; }

// ----------
// Property
// ----------

inline const std::string& Property::Name() const
{ return Attribute("Name"); }

inline std::string& Property::Name()
{ return Attribute("Name"); }

inline Property& Property::Name(const std::string& name)
{ Attribute("Name", name); return *this; }

inline const std::string& Property::Operator() const
{ return Attribute("Operator"); }

inline std::string& Property::Operator()
{ return Attribute("Operator"); }

inline Property& Property::Operator(const std::string& op)
{ Attribute("Operator", op); return *this; }

inline const std::string& Property::Value() const
{ return Attribute("Value"); }

inline std::string& Property::Value()
{ return Attribute("Value"); }

inline Property& Property::Value(const std::string& value)
{ Attribute("Value", value); return *this; }

// ------------
// Provenance
// ------------

inline const std::string& Provenance::CreatedBy() const
{ return Attribute("CreatedBy"); }

inline std::string& Provenance::CreatedBy()
{ return Attribute("CreatedBy"); }

inline Provenance& Provenance::CreatedBy(const std::string& createdBy)
{ Attribute("CreatedBy", createdBy); return *this; }

inline const std::string& Provenance::CommonServicesInstanceId() const
{ return ChildText("CommonServicesInstanceId"); }

inline std::string& Provenance::CommonServicesInstanceId()
{ return ChildText("CommonServicesInstanceId"); }

inline Provenance& Provenance::CommonServicesInstanceId(const std::string& id)
{ ChildText("CommonServicesInstanceId", id); return *this; }

inline const std::string& Provenance::CreatorUserId() const
{ return ChildText("CreatorUserId"); }

inline std::string& Provenance::CreatorUserId()
{ return ChildText("CreatorUserId"); }

inline Provenance& Provenance::CreatorUserId(const std::string& id)
{ ChildText("CreatorUserId", id); return *this; }

inline const std::string& Provenance::ParentJobId() const
{ return ChildText("ParentJobId"); }

inline std::string& Provenance::ParentJobId()
{ return ChildText("ParentJobId"); }

inline Provenance& Provenance::ParentJobId(const std::string& id)
{ ChildText("ParentJobId", id); return *this; }

inline Provenance& Provenance::ParentTool(const PacBio::BAM::ParentTool& tool)
{ ParentTool() = tool; return *this; }

} // namespace BAM
} // namespace PacBio
