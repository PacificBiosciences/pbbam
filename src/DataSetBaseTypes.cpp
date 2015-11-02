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

#include "pbbam/DataSetTypes.h"
#include "pbbam/internal/DataSetBaseTypes.h"
#include "DataSetUtils.h"
#include "TimeUtils.h"
#include <boost/algorithm/string.hpp>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace std;

// ----------------
// BaseEntityType
// ----------------

BaseEntityType::BaseEntityType(const std::string& label, const XsdType& xsd)
    : DataSetElement(label, xsd)
{
    if (Version().empty())
        Version(internal::XML_VERSION);
}

DEFINE_ACCESSORS(BaseEntityType, Extensions, Extensions)

BaseEntityType& BaseEntityType::Extensions(const PacBio::BAM::Extensions& extensions)
{ Extensions() = extensions; return *this; }

// ----------------
// DataEntityType
// ----------------

DataEntityType::DataEntityType(const std::string& label, const XsdType& xsd)
    : BaseEntityType(label, xsd)
{ }

// -----------------
// IndexedDataType
// -----------------

IndexedDataType::IndexedDataType(const string& metatype,
                                 const string& filename,
                                 const string& label, 
                                 const XsdType &xsd)
    : InputOutputDataType(metatype, filename, label, xsd)
{ }

DEFINE_ACCESSORS(IndexedDataType, FileIndices, FileIndices)

IndexedDataType& IndexedDataType::FileIndices(const PacBio::BAM::FileIndices& indices)
{ FileIndices() = indices; return *this; }

// ---------------------
// InputOutputDataType
// ---------------------

InputOutputDataType::InputOutputDataType(const string& metatype,
                                         const string& filename,
                                         const string& label,
                                         const XsdType &xsd)
    : StrictEntityType(metatype, label, xsd)
{  
    ResourceId(filename);
}

// ----------------
// StrictEntityType
// ----------------

StrictEntityType::StrictEntityType(const string& metatype, 
                                   const string& label, 
                                   const XsdType& xsd)
    : BaseEntityType(label, xsd)
{ 
    // MetaType
    MetaType(metatype);

    // TimeStampedName
    const size_t numChars = metatype.size();
    string transformedMetatype;
    transformedMetatype.resize(numChars);
    for (size_t i = 0; i < numChars; ++i) {
        const char c = metatype.at(i);
        transformedMetatype[i] = ((c == '.') ? '_' : tolower(c));
    }
    const string& tsn = transformedMetatype + "-" + internal::ToDataSetFormat(internal::CurrentTime());
    TimeStampedName(tsn);

    // UniqueId
    UniqueId(internal::GenerateUuid());
}
