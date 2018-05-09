// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "DataSetIO.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <exception>
#include <fstream>
#include <iostream>

#include <boost/algorithm/string.hpp>

#include "FileUtils.h"
#include "FofnReader.h"
#include "StringUtils.h"
#include "XmlReader.h"
#include "XmlWriter.h"
#include "pbbam/MakeUnique.h"

namespace PacBio {
namespace BAM {

using DataSetPtr = std::shared_ptr<DataSetBase>;

namespace internal {

static std::unique_ptr<DataSetBase> FromXml(const std::string& xmlFn)
{
    std::ifstream in(xmlFn);
    if (!in) throw std::runtime_error{"could not open XML file for reading: " + xmlFn};
    return XmlReader::FromStream(in);
}

static std::unique_ptr<DataSetBase> FromBam(const std::string& bamFn)
{
    // peek at sort order to determine if file should be an AlignmentSet or else SubreadSet
    const auto bamFile = BamFile{bamFn};
    const auto& header = bamFile.Header();
    const auto aligned = header.SortOrder() == "coordinate";

    std::unique_ptr<DataSetBase> dataset;
    if (aligned)
        dataset = std::make_unique<AlignmentSet>();
    else
        dataset = std::make_unique<SubreadSet>();

    auto& resources = dataset->ExternalResources();
    resources.Add(ExternalResource(BamFile(bamFn)));
    return dataset;
}

static std::unique_ptr<DataSetBase> FromFasta(const std::string& fasta)
{
    // make FASTA data set
    auto dataset = std::make_unique<ReferenceSet>();
    auto& resources = dataset->ExternalResources();
    resources.Add(ExternalResource("PacBio.ReferenceFile.ReferenceFastaFile", fasta));
    return std::move(dataset);
}

static std::unique_ptr<DataSetBase> FromFofn(const std::string& fofn)
{
    const auto fofnDir = FileUtils::DirectoryName(fofn);
    std::ifstream in(fofn);
    if (!in) throw std::runtime_error{"could not open FOFN for reading: " + fofn};

    auto filenames = FofnReader::Files(in);
    std::transform(
        filenames.begin(), filenames.end(), filenames.begin(),
        [&fofnDir](const std::string fn) { return FileUtils::ResolvedFilePath(fn, fofnDir); });
    return DataSetIO::FromUris(filenames);
}

static std::unique_ptr<DataSetBase> FromUri(const std::string& uri)
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
    else if (boost::algorithm::iends_with(uri, ".fasta") ||
             boost::algorithm::iends_with(uri, ".fa")) {
        return FromFasta(uri);
    }

    // unknown filename extension
    throw std::runtime_error{"unsupported extension on input file: " + uri};
}

std::unique_ptr<DataSetBase> DataSetIO::FromUri(const std::string& uri)
{
    return FromUris(std::vector<std::string>(1, uri));
}

std::unique_ptr<DataSetBase> DataSetIO::FromUris(const std::vector<std::string>& uris)
{
    if (uris.empty())
        throw std::runtime_error{"empty input URI list"};  // or just return empty, generic DataSet?

    // create dataset(s) from URI(s)
    std::vector<std::unique_ptr<DataSetBase> > datasets;
    datasets.reserve(uris.size());
    for (const auto& uri : uris)
        datasets.emplace_back(internal::FromUri(uri));
    assert(!datasets.empty());

    // if only 1, just return
    if (datasets.size() == 1) return std::unique_ptr<DataSetBase>(datasets.front().release());

    // else merge
    else {
        auto& result = datasets.front();
        for (const auto& dataset : datasets)
            *result += *dataset;
        return std::move(result);
    }
}

std::unique_ptr<DataSetBase> DataSetIO::FromXmlString(const std::string& xml)
{
    if (xml.empty()) throw std::runtime_error{"empty XML string"};
    std::istringstream s{xml};
    return XmlReader::FromStream(s);
}

void DataSetIO::ToFile(const std::unique_ptr<DataSetBase>& dataset, const std::string& fn)
{
    std::ofstream out(fn);
    if (!out) throw std::runtime_error{"could not open XML file for writing: " + fn};
    XmlWriter::ToStream(dataset, out);
}

void DataSetIO::ToStream(const std::unique_ptr<DataSetBase>& dataset, std::ostream& out)
{
    XmlWriter::ToStream(dataset, out);
}

}  // namespace internal
}  // namespace BAM
}  // namespace PacBio
