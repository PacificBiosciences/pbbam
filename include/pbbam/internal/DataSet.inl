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
/// \file DataSet.inl
/// \brief Inline implementations for the DataSet class.
//
// Author: Derek Barnett

#include "pbbam/DataSet.h"

namespace PacBio {
namespace BAM {

inline const std::string& DataSet::Attribute(const std::string& name) const
{ return d_->Attribute(name); }

inline std::string& DataSet::Attribute(const std::string& name)
{ return d_->Attribute(name); }

inline DataSet& DataSet::Attribute(const std::string& name, const std::string& value)
{ d_->Attribute(name, value); return *this; }

inline const std::string& DataSet::CreatedAt(void) const
{ return d_->CreatedAt(); }

inline std::string& DataSet::CreatedAt(void)
{ return d_->CreatedAt(); }

inline DataSet& DataSet::CreatedAt(const std::string& createdAt)
{ d_->CreatedAt(createdAt); return *this; }

inline const PacBio::BAM::Extensions& DataSet::Extensions(void) const
{ return d_->Extensions(); }

inline PacBio::BAM::Extensions& DataSet::Extensions(void)
{ return d_->Extensions(); }

inline DataSet& DataSet::Extensions(const PacBio::BAM::Extensions& extensions)
{ d_->Extensions(extensions); return *this; }

inline const PacBio::BAM::ExternalResources& DataSet::ExternalResources(void) const
{ return d_->ExternalResources(); }

inline PacBio::BAM::ExternalResources& DataSet::ExternalResources(void)
{ return d_->ExternalResources(); }

inline DataSet& DataSet::ExternalResources(const PacBio::BAM::ExternalResources& resources)
{ d_->ExternalResources(resources); return *this; }

inline const PacBio::BAM::Filters& DataSet::Filters(void) const
{ return d_->Filters(); }

inline PacBio::BAM::Filters& DataSet::Filters(void)
{ return d_->Filters(); }

inline DataSet& DataSet::Filters(const PacBio::BAM::Filters& filters)
{ d_->Filters(filters); return *this; }

inline const std::string& DataSet::Format(void) const
{ return d_->Format(); }

inline std::string& DataSet::Format(void)
{ return d_->Format(); }

inline DataSet& DataSet::Format(const std::string& format)
{ d_->Format(format); return *this; }

inline const PacBio::BAM::DataSetMetadata& DataSet::Metadata(void) const
{ return d_->Metadata(); }

inline PacBio::BAM::DataSetMetadata& DataSet::Metadata(void)
{ return d_->Metadata(); }

inline DataSet& DataSet::Metadata(const PacBio::BAM::DataSetMetadata& metadata)
{ d_->Metadata(metadata); return *this; }

inline const std::string& DataSet::MetaType(void) const
{ return d_->MetaType(); }

inline std::string& DataSet::MetaType(void)
{ return d_->MetaType(); }

inline DataSet& DataSet::MetaType(const std::string& metatype)
{ d_->MetaType(metatype); return *this; }

inline const std::string& DataSet::ModifiedAt(void) const
{ return d_->ModifiedAt(); }

inline std::string& DataSet::ModifiedAt(void)
{ return d_->ModifiedAt(); }

inline DataSet& DataSet::ModifiedAt(const std::string& modifiedAt)
{ d_->ModifiedAt(modifiedAt); return *this; }

inline const std::string& DataSet::Name(void) const
{ return d_->Name(); }

inline std::string& DataSet::Name(void)
{ return d_->Name(); }

inline DataSet& DataSet::Name(const std::string& name)
{ d_->Name(name); return *this; }

inline const std::string& DataSet::ResourceId(void) const
{ return d_->ResourceId(); }

inline std::string& DataSet::ResourceId(void)
{ return d_->ResourceId(); }

inline DataSet& DataSet::ResourceId(const std::string& resourceId)
{ d_->ResourceId(resourceId); return *this; }

inline const PacBio::BAM::SubDataSets& DataSet::SubDataSets(void) const
{ return d_->SubDataSets(); }

inline PacBio::BAM::SubDataSets& DataSet::SubDataSets(void)
{ return d_->SubDataSets(); }

inline DataSet& DataSet::SubDataSets(const PacBio::BAM::SubDataSets& subdatasets)
{ d_->SubDataSets(subdatasets); return *this; }

inline const std::string& DataSet::Tags(void) const
{ return d_->Tags(); }

inline std::string& DataSet::Tags(void)
{ return d_->Tags(); }

inline DataSet& DataSet::Tags(const std::string& tags)
{ d_->Tags(tags); return *this; }

inline const std::string& DataSet::TimeStampedName(void) const
{ return d_->TimeStampedName(); }

inline std::string& DataSet::TimeStampedName(void)
{ return d_->TimeStampedName(); }

inline DataSet& DataSet::TimeStampedName(const std::string& timeStampedName)
{ d_->TimeStampedName(timeStampedName); return *this; }

inline PacBio::BAM::DataSet::TypeEnum DataSet::Type(void) const
{ return DataSet::NameToType(TypeName()); }

inline DataSet& DataSet::Type(const DataSet::TypeEnum type)
{ d_->Label(DataSet::TypeToName(type)); return *this; }

inline std::string DataSet::TypeName(void) const
{ return d_->LocalNameLabel().to_string(); }

inline const std::string& DataSet::UniqueId(void) const
{ return d_->UniqueId(); }

inline std::string& DataSet::UniqueId(void)
{ return d_->UniqueId(); }

inline DataSet& DataSet::UniqueId(const std::string& uuid)
{ d_->UniqueId(uuid); return *this; }

inline const std::string& DataSet::Version(void) const
{ return d_->Version(); }

inline std::string& DataSet::Version(void)
{ return d_->Version(); }

inline DataSet& DataSet::Version(const std::string& version)
{ d_->Version(version); return *this; }

} // namespace BAM
} // namespace PacBio
