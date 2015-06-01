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

#include "pbbam/dataset/DataSet.h"
#include "pbbam/dataset/DataSetMetadata.h"
#include "DataSetIO.h"
#include "StringUtils.h"
#include <boost/algorithm/string.hpp>
#include <set>
#include <vector>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

// empty, "null" components

static
const DataSetMetadata& NullMetadata(void)
{
    static const DataSetMetadata empty;
    return empty;
}

} // namespace internal
} // namespace BAM
} // namespace PacBio

DataSet::DataSet(void)
    : DataSetBase(DataSetType::GENERIC)
{ }

DataSet::DataSet(const BamFile& bamFile)
    : DataSetBase(DataSetType::GENERIC)
{ DataSetBase::operator=(std::move(DataSetIO::FromUri(bamFile.Filename()))); }

DataSet& DataSet::operator+=(const DataSet& other)
{
    // fail on conflicting metadata, just for now to simplify
    const DataSetMetadata& metadata = Metadata();
    const DataSetMetadata& otherMetadata = other.Metadata();
    if (metadata != otherMetadata)
        throw std::exception();
    DataSetBase::operator+=(other);
    return *this;
}

vector<BamFile> DataSet::BamFiles(void) const
{ return ExternalDataReferenceList().BamFiles(); }

const DataSetMetadata& DataSet::Metadata(void) const
{
    try {
        return Child<DataSetMetadata>("DataSetMetadata");
    } catch (std::exception&) {
        return internal::NullMetadata();
    }
}

DataSetMetadata& DataSet::Metadata(void)
{
    if (!HasChild("DataSetMetadata"))
        AddChild(internal::NullMetadata());
    return Child<DataSetMetadata>("DataSetMetadata");
}

AlignmentSet DataSet::ToAlignmentSet(void) const
{ return AlignmentSet(*this); }

BarcodeSet DataSet::ToBarcodeSet(void) const
{ return BarcodeSet(*this); }

CcsReadSet DataSet::ToCcsReadSet(void) const
{ return CcsReadSet(*this); }

ContigSet DataSet::ToContigSet(void) const
{ return ContigSet(*this); }

ReferenceSet DataSet::ToReferenceSet(void) const
{  return ReferenceSet(*this); }

SubreadSet DataSet::ToSubreadSet(void) const
{ return SubreadSet(*this); }
