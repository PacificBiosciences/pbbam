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

#include "DataSetIO.h"
#include "FofnReader.h"
#include "StringUtils.h"
#include "XmlReader.h"
#include "XmlWriter.h"
#include <boost/algorithm/string.hpp>
#include <exception>
#include <fstream>
#include <iostream>
#include <cassert>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

static
DataSetBase FromXml(const string& xmlFn)
{
    ifstream in(xmlFn);
    if (!in)
        throw std::exception();
    return XmlReader::FromStream(in);
}

static
DataSetBase FromBam(const string& bamFn)
{
    DataSet dataset(DataSetType::SUBREADSET);
    BamFile bamFile(bamFn);
    dataset.AddExternalDataReference(ExternalDataReference(bamFile));
    return dataset;
}

static
DataSetBase FromFofn(const string& fofn)
{
    ifstream in(fofn);
    if (!in)
        throw std::exception();
    const vector<string> filenames = std::move(FofnReader::Files(in));
    return DataSetIO::FromUris(filenames);
}

static
DataSetBase FromUri(const string& uri)
{
    // NOTE: this says URI, but we're not quite handling filenames as true URIs
    //       basically just treating as a regular filename for now

    // handle on extension
    if (boost::algorithm::iends_with(uri, ".xml"))
        return FromXml(uri);
    else if (boost::algorithm::iends_with(uri, ".bam"))
        return FromBam(uri);
    else if (boost::algorithm::iends_with(uri, ".fofn"))
        return FromFofn(uri);

    // unknown filename extension
    throw std::exception();
}

} // namespace internal
} // namespace BAM
} // namespace PacBio

DataSetBase DataSetIO::FromUri(const string& uri)
{ return FromUris(vector<string>(1, uri)); }

DataSetBase DataSetIO::FromUris(const vector<string>& uris)
{
    if (uris.empty())
        return DataSetBase(); // exception? or is quiet, empty dataset OK?

    vector<DataSetBase> datasets;
    datasets.reserve(uris.size());
    for ( const auto& uri : uris )
        datasets.push_back(internal::FromUri(uri));
    assert(!datasets.empty());

    // merge datasets (if more than 1 loaded)
    DataSetBase result = datasets.front();
    for ( size_t i = 1; i < datasets.size(); ++i)
        result += datasets.at(i);
    return result;
}

void DataSetIO::ToFile(const DataSetBase& dataset, const string& fn)
{
    ofstream out(fn);
    if (!out)
        throw std::exception();
    XmlWriter::ToStream(dataset, out);
}

void DataSetIO::ToStream(const DataSetBase& dataset, ostream &out)
{ XmlWriter::ToStream(dataset, out); }

