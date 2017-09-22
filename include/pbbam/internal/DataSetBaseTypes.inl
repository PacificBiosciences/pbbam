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

#include "pbbam/internal/DataSetBaseTypes.h"

namespace PacBio {
namespace BAM {
namespace internal {

// ----------------
// BaseEntityType
// ----------------

inline const std::string& BaseEntityType::CreatedAt() const
{ return Attribute("CreatedAt"); }

inline std::string& BaseEntityType::CreatedAt()
{ return Attribute("CreatedAt"); }

inline BaseEntityType& BaseEntityType::CreatedAt(const std::string& createdAt)
{ Attribute("CreatedAt", createdAt); return *this; }

inline const std::string& BaseEntityType::Description() const
{ return Attribute("Description"); }

inline std::string& BaseEntityType::Description()
{ return Attribute("Description"); }

inline BaseEntityType& BaseEntityType::Description(const std::string& description)
{ Attribute("Description", description); return *this; }

inline const std::string& BaseEntityType::Format() const
{ return Attribute("Format"); }

inline std::string& BaseEntityType::Format()
{ return Attribute("Format"); }

inline BaseEntityType& BaseEntityType::Format(const std::string& format)
{ Attribute("Format", format); return *this; }

inline const std::string& BaseEntityType::ModifiedAt() const
{ return Attribute("ModifiedAt"); }

inline std::string& BaseEntityType::ModifiedAt()
{ return Attribute("ModifiedAt"); }

inline BaseEntityType& BaseEntityType::ModifiedAt(const std::string& modifiedAt)
{ Attribute("ModifiedAt", modifiedAt); return *this; }

inline const std::string& BaseEntityType::Name() const
{ return Attribute("Name"); }

inline std::string& BaseEntityType::Name()
{ return Attribute("Name"); }

inline BaseEntityType& BaseEntityType::Name(const std::string& name)
{ Attribute("Name", name); return *this; }

inline const std::string& BaseEntityType::ResourceId() const
{ return Attribute("ResourceId"); }

inline std::string& BaseEntityType::ResourceId()
{ return Attribute("ResourceId"); }

inline BaseEntityType& BaseEntityType::ResourceId(const std::string& resourceId)
{ Attribute("ResourceId", resourceId); return *this; }

inline const std::string& BaseEntityType::Tags() const
{ return Attribute("Tags"); }

inline std::string& BaseEntityType::Tags()
{ return Attribute("Tags"); }

inline BaseEntityType& BaseEntityType::Tags(const std::string& tags)
{ Attribute("Tags", tags); return *this; }

inline const std::string& BaseEntityType::Version() const
{ return Attribute("Version"); }

inline std::string& BaseEntityType::Version()
{ return Attribute("Version"); }

inline BaseEntityType& BaseEntityType::Version(const std::string& version)
{ Attribute("Version", version); return *this; }

// ----------------
// DataEntityType
// ----------------

inline const std::string& DataEntityType::Checksum() const
{ return ChildText("Checksum"); }

inline std::string& DataEntityType::Checksum()
{ return ChildText("Checksum"); }

inline DataEntityType& DataEntityType::Checksum(const std::string& checksum)
{ ChildText("Checksum", checksum); return *this; }

inline const std::string& DataEntityType::EncodedValue() const
{ return ChildText("EncodedValue"); }

inline std::string& DataEntityType::EncodedValue()
{ return ChildText("EncodedValue"); }

inline DataEntityType& DataEntityType::EncodedValue(const std::string& encodedValue)
{ ChildText("EncodedValue", encodedValue); return *this; }

inline const std::string& DataEntityType::MetaType() const
{ return Attribute("MetaType"); }

inline std::string& DataEntityType::MetaType()
{ return Attribute("MetaType"); }

inline DataEntityType& DataEntityType::MetaType(const std::string& metatype)
{ Attribute("MetaType", metatype); return *this; }

inline const std::string& DataEntityType::SimpleValue() const
{ return Attribute("SimpleValue"); }

inline std::string& DataEntityType::SimpleValue()
{ return Attribute("SimpleValue"); }

inline DataEntityType& DataEntityType::SimpleValue(const std::string& simpleValue)
{ Attribute("SimpleValue", simpleValue); return *this; }

inline const std::string& DataEntityType::TimeStampedName() const
{ return Attribute("TimeStampedName"); }

inline std::string& DataEntityType::TimeStampedName()
{ return Attribute("TimeStampedName"); }

inline DataEntityType& DataEntityType::TimeStampedName(const std::string& timeStampedName)
{ Attribute("TimeStampedName", timeStampedName); return *this; }

inline const std::string& DataEntityType::UniqueId() const
{ return Attribute("UniqueId"); }

inline std::string& DataEntityType::UniqueId()
{ return Attribute("UniqueId"); }

inline DataEntityType& DataEntityType::UniqueId(const std::string& uuid)
{ Attribute("UniqueId", uuid); return *this; }

inline const std::string& DataEntityType::ValueDataType() const
{ return Attribute("ValueDataType"); }

inline std::string& DataEntityType::ValueDataType()
{ return Attribute("ValueDataType"); }

inline DataEntityType& DataEntityType::ValueDataType(const std::string& valueDataType)
{ Attribute("ValueDataType", valueDataType); return *this; }

// ----------------
// StrictEntityType
// ----------------

inline const std::string& StrictEntityType::MetaType() const
{ return Attribute("MetaType"); }

inline std::string& StrictEntityType::MetaType()
{ return Attribute("MetaType"); }

inline StrictEntityType& StrictEntityType::MetaType(const std::string& metatype)
{ Attribute("MetaType", metatype); return *this; }

inline const std::string& StrictEntityType::TimeStampedName() const
{ return Attribute("TimeStampedName"); }

inline std::string& StrictEntityType::TimeStampedName()
{ return Attribute("TimeStampedName"); }

inline StrictEntityType& StrictEntityType::TimeStampedName(const std::string& timeStampedName)
{ Attribute("TimeStampedName", timeStampedName); return *this; }

inline const std::string& StrictEntityType::UniqueId() const
{ return Attribute("UniqueId"); }

inline std::string& StrictEntityType::UniqueId()
{ return Attribute("UniqueId"); }

inline StrictEntityType& StrictEntityType::UniqueId(const std::string& uuid)
{ Attribute("UniqueId", uuid); return *this; }

} // namespace internal
} // namespace BAM
} // namespace PacBio
