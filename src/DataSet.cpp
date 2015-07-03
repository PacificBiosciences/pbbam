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

#include "pbbam/DataSet.h"
#include "pbbam/DataSetTypes.h"
#include "pbbam/internal/DataSetBaseTypes.h"
#include "DataSetIO.h"
#include <unordered_map>
#include <iostream>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace std;

DataSet::DataSet(void)
    : d_(new DataSetBase)
{ }

DataSet::DataSet(const DataSet::TypeEnum type)
    : d_(nullptr)
{
    switch(type) {
        case DataSet::GENERIC             : d_.reset(new DataSetBase); break;
        case DataSet::ALIGNMENT           : d_.reset(new AlignmentSet); break;
        case DataSet::BARCODE             : d_.reset(new BarcodeSet); break;
        case DataSet::CONSENSUS_ALIGNMENT : d_.reset(new ConsensusAlignmentSet); break;
        case DataSet::CONSENSUS_READ      : d_.reset(new ConsensusReadSet); break;
        case DataSet::CONTIG              : d_.reset(new ContigSet); break;
        case DataSet::HDF_SUBREAD         : d_.reset(new HdfSubreadSet); break;
        case DataSet::REFERENCE           : d_.reset(new ReferenceSet); break;
        case DataSet::SUBREAD             : d_.reset(new SubreadSet); break;
        default:
            throw std::runtime_error("unsupported dataset type"); // unknown type
    }
}

DataSet::DataSet(const BamFile& bamFile)
    : d_(internal::DataSetIO::FromUri(bamFile.Filename()))
{ }

DataSet::DataSet(const string& filename)
    : d_(internal::DataSetIO::FromUri(filename))
{ }

DataSet::DataSet(const DataSet& other)
{
    DataSetBase* otherDataset = other.d_.get();
    DataSetElement* copyDataset = new DataSetElement(*otherDataset);
    d_.reset(static_cast<DataSetBase*>(copyDataset));
}

DataSet::DataSet(DataSet&& other)
    : d_(std::move(other.d_))
{
    assert(other.d_.get() == nullptr);
}

DataSet& DataSet::operator=(const DataSet& other)
{
    DataSetBase* otherDataset = other.d_.get();
    DataSetElement* copyDataset = new DataSetElement(*otherDataset);
    d_.reset(static_cast<DataSetBase*>(copyDataset));
    return *this;
}

DataSet& DataSet::operator=(DataSet&& other)
{
    d_ = std::move(other.d_);
    return *this;
}

DataSet::~DataSet(void) { }

DataSet& DataSet::operator+=(const DataSet&)
{
    return *this;
}

bool DataSet::operator==(const DataSet&) const
{
    return true;
}

DataSet::TypeEnum DataSet::NameToType(const string& typeName)
{
    static std::unordered_map<std::string, DataSet::TypeEnum> lookup;
    if (lookup.empty()) {
        lookup["DataSet"] = DataSet::GENERIC;
        lookup["AlignmentSet"] = DataSet::ALIGNMENT;
        lookup["BarcodeSet"] = DataSet::BARCODE;
        lookup["ConsensusAlignmentSet"] = DataSet::CONSENSUS_ALIGNMENT;
        lookup["ConsensusReadSet"] = DataSet::CONSENSUS_READ;
        lookup["ContigSet"] = DataSet::CONTIG;
        lookup["HdfSubreadSet"] = DataSet::HDF_SUBREAD;
        lookup["ReferenceSet"] = DataSet::REFERENCE;
        lookup["SubreadSet"] = DataSet::SUBREAD;
    }
    return lookup.at(typeName); // throws if unknown typename
}

void DataSet::Save(const std::string& outputFilename)
{ internal::DataSetIO::ToFile(d_, outputFilename); }

void DataSet::SaveToStream(ostream& out)
{ internal::DataSetIO::ToStream(d_, out); }

string DataSet::TypeToName(const DataSet::TypeEnum& type)
{
    switch(type) {
        case DataSet::GENERIC             : return "DataSet";
        case DataSet::ALIGNMENT           : return "AlignmentSet";
        case DataSet::BARCODE             : return "BarcodeSet";
        case DataSet::CONSENSUS_ALIGNMENT : return "ConsensusAlignmentSet";
        case DataSet::CONSENSUS_READ      : return "ConsensusReadSet";
        case DataSet::CONTIG              : return "ContigSet";
        case DataSet::HDF_SUBREAD         : return "HdfSubreadSet";
        case DataSet::REFERENCE           : return "ReferenceSet";
        case DataSet::SUBREAD             : return "SubreadSet";
        default:
            throw std::runtime_error("unsupported dataset type"); // unknown type
    }
}

//DataSet& DataSet::operator+=(const DataSet& other)
//{
//    // fail on conflicting metadata, just for now to simplify
//    const DataSetMetadata& metadata = Metadata();
//    const DataSetMetadata& otherMetadata = other.Metadata();
//    if (metadata != otherMetadata)
//        throw std::exception();
//    DataSetBase::operator+=(other);
//    return *this;
//}

