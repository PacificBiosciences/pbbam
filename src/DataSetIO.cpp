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
#include "FileUtils.h"
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

typedef std::shared_ptr<DataSetBase> DataSetPtr;

namespace PacBio {
namespace BAM {
namespace internal {

static
unique_ptr<DataSetBase> FromXml(const string& xmlFn)
{
    ifstream in(xmlFn);
    if (!in)
        throw std::runtime_error("could not open XML file for reading");
    return XmlReader::FromStream(in);
}

static
unique_ptr<DataSetBase> FromBam(const string& bamFn)
{
    unique_ptr<DataSetBase> dataset(new SubreadSet);
    ExternalResources& resources = dataset->ExternalResources();
    resources.Add(ExternalResource(BamFile(bamFn)));
    return dataset;
}

static
unique_ptr<DataSetBase> FromFofn(const string& fofn)
{
    const string fofnDir = internal::FileUtils::DirectoryName(fofn);
    ifstream in(fofn);
    if (!in)
        throw std::runtime_error("could not open FOFN for reading");

    vector<string> filenames = std::move(FofnReader::Files(in));
    for (size_t i = 0; i < filenames.size(); ++i)
        filenames[i] = internal::FileUtils::ResolvedFilePath(filenames[i], fofnDir);
    return DataSetIO::FromUris(filenames);
}

static
unique_ptr<DataSetBase> FromUri(const string& uri)
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
    throw std::runtime_error("unsupported input file extension");
}

} // namespace internal
} // namespace BAM
} // namespace PacBio

std::unique_ptr<DataSetBase> DataSetIO::FromUri(const std::string& uri)
{
    return FromUris(vector<string>(1, uri));
}

std::unique_ptr<DataSetBase> DataSetIO::FromUris(const std::vector<std::string>& uris)
{
    if (uris.empty())
        throw std::runtime_error("empty input URI list"); // or just return empty, generic DataSet?

    // create dataset(s) from URI(s)
    vector< unique_ptr<DataSetBase> > datasets;
    datasets.reserve(uris.size());
    for ( const auto& uri : uris )
        datasets.push_back(internal::FromUri(uri));
    assert(!datasets.empty());

    // if only 1, just return
    if (datasets.size() == 1)
        return std::unique_ptr<DataSetBase>(datasets.front().release());

    // else merge
    else {
        unique_ptr<DataSetBase>& result = datasets.front();
        for (size_t i = 1; i < datasets.size(); ++i)
            *result += *datasets.at(i);
        return unique_ptr<DataSetBase>(result.release());
    }
}

std::unique_ptr<DataSetBase> DataSetIO::FromXmlString(const string& xml)
{
    if (xml.empty())
        throw std::runtime_error("empty XML string");
    stringstream s(xml);
    return XmlReader::FromStream(s);
}

void DataSetIO::ToFile(const std::unique_ptr<DataSetBase>& dataset,
                       const string& fn)
{
    ofstream out(fn);
    if (!out)
        throw std::runtime_error("could not open XML for writing");
    XmlWriter::ToStream(dataset, out);
}

void DataSetIO::ToStream(const std::unique_ptr<DataSetBase>& dataset, ostream &out)
{ XmlWriter::ToStream(dataset, out); }
