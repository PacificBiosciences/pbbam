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
using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace std;

// ------------------
// AlignmentSetType
// ------------------

//AlignmentSetType::AlignmentSetType(const std::string& label)
//    : DataSetType(label)
//{ }

// ----------------
// BarcodeSetType
// ----------------

//BarcodeSetType::BarcodeSetType(const std::string& label)
//    : DataSetType(label)
//{ }

// ----------------
// BaseEntityType
// ----------------

BaseEntityType::BaseEntityType(const std::string& label)
    : DataSetElement(label)
{ }

DEFINE_ACCESSORS(BaseEntityType, Extensions, Extensions)

BaseEntityType& BaseEntityType::Extensions(const PacBio::BAM::Extensions& extensions)
{ Extensions() = extensions; return *this; }

//ContigSetType::ContigSetType(const std::string& label)
//    : DataSetType(label)
//{ }

// ----------------
// DataEntityType
// ----------------

DataEntityType::DataEntityType(const std::string& label)
    : BaseEntityType(label)
{ }

//// ---------------------
//// DataSetMetadataType
//// ---------------------

//DataSetMetadataType::DataSetMetadataType(const std::string& numRecords,
//                                         const std::string& totalLength)
//    : DataSetElement("DataSetMetadata")
//{
//    NumRecords(numRecords);
//    TotalLength(totalLength);
//}

//DEFINE_ACCESSORS(DataSetMetadataType, Provenance, Provenance)

//DataSetMetadataType& DataSetMetadataType::Provenance(const PacBio::BAM::Provenance& provenance)
//{ Provenance() = provenance; return *this; }

//// -------------
//// DataSetType
//// -------------

//DataSetType::DataSetType(const std::string& label)
//    : BaseEntityType(label)
//{ }

//DEFINE_ACCESSORS(DataSetType, ExternalResources, ExternalResources)
//DEFINE_ACCESSORS(DataSetType, Filters,           Filters)
//DEFINE_ACCESSORS(DataSetType, DataSetMetadata,   Metadata)

//DataSetType& DataSetType::ExternalResources(const PacBio::BAM::ExternalResources& resources)
//{ ExternalResources() = resources; return *this; }

//DataSetType& DataSetType::Filters(const PacBio::BAM::Filters& filters)
//{ Filters() = filters; return *this; }

//DataSetType& DataSetType::Metadata(const DataSetMetadata& metadata)
//{ Metadata() = metadata; return *this; }


// -----------------
// IndexedDataType
// -----------------

IndexedDataType::IndexedDataType(const std::string& label)
    : InputOutputDataType(label)
{ }

DEFINE_ACCESSORS(IndexedDataType, FileIndices, FileIndices)

IndexedDataType& IndexedDataType::FileIndices(const PacBio::BAM::FileIndices& indices)
{ FileIndices() = indices; return *this; }

// ---------------------
// InputOutputDataType
// ---------------------

InputOutputDataType::InputOutputDataType(const std::string& label)
    : DataEntityType(label)
{ }

// -------------
// ReadSetType
// -------------

//ReadSetType::ReadSetType(const std::string& label)
//    : DataSetType(label)
//{ }

//// ----------------
//// SubreadSetType
//// ----------------

//SubreadSetType::SubreadSetType(const std::string& label)
//    : ReadSetType(label)
//{ }
