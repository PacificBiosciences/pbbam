// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "DataSetIO.h"

#include <cassert>
#include <cstddef>

#include <algorithm>
#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <boost/algorithm/string.hpp>

#include "pbbam/StringUtilities.h"

#include "FileUtils.h"
#include "FofnReader.h"
#include "XmlReader.h"
#include "XmlWriter.h"

namespace PacBio {
namespace BAM {
namespace {

struct DataSetFileException : public std::exception
{
    DataSetFileException(std::string filename, std::string reason) : std::exception{}
    {
        std::ostringstream s;
        s << "[pbbam] dataset I/O ERROR: " << reason << ":\n"
          << "  file: " << filename;
        msg_ = s.str();
    }

    const char* what() const noexcept override { return msg_.c_str(); }

    std::string msg_;
};

std::unique_ptr<DataSetBase> DataSetFromXml(const std::string& xmlFn)
{
    std::ifstream in(xmlFn);
    if (!in) throw DataSetFileException{xmlFn, "could not open XML file for reading"};
    return XmlReader::FromStream(in);
}

std::unique_ptr<DataSetBase> DataSetFromBam(const std::string& bamFn)
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

std::unique_ptr<DataSetBase> DataSetFromFasta(const std::string& fasta)
{
    // make FASTA data set
    auto dataset = std::make_unique<ReferenceSet>();
    auto& resources = dataset->ExternalResources();
    resources.Add(ExternalResource("PacBio.ReferenceFile.ReferenceFastaFile", fasta));
    return
#ifdef __INTEL_COMPILER
        std::move(
#endif
            dataset
#ifdef __INTEL_COMPILER
            )
#endif
            ;
}

std::unique_ptr<DataSetBase> DataSetFromFofn(const std::string& fofn)
{
    const auto fofnDir = FileUtils::DirectoryName(fofn);
    std::ifstream in(fofn);
    if (!in) throw DataSetFileException{fofn, "could not open FOFN for reading"};

    auto filenames = FofnReader::Files(in);
    std::transform(
        filenames.begin(), filenames.end(), filenames.begin(),
        [&fofnDir](const std::string fn) { return FileUtils::ResolvedFilePath(fn, fofnDir); });
    return DataSetIO::FromUris(filenames);
}

std::unique_ptr<DataSetBase> DataSetFromUri(const std::string& uri)
{
    // NOTE: this says URI, but we're not quite handling filenames as true URIs
    //       basically just treating as a regular filename for now

    // handle on extension
    if (boost::algorithm::iends_with(uri, ".xml"))
        return DataSetFromXml(uri);
    else if (boost::algorithm::iends_with(uri, ".bam"))
        return DataSetFromBam(uri);
    else if (boost::algorithm::iends_with(uri, ".fofn"))
        return DataSetFromFofn(uri);
    else if (boost::algorithm::iends_with(uri, ".fasta") ||
             boost::algorithm::iends_with(uri, ".fa")) {
        return DataSetFromFasta(uri);
    }

    // unknown filename extension
    throw DataSetFileException{uri, "unsupported extension on input file"};
}

}  // namespace

std::unique_ptr<DataSetBase> DataSetIO::FromUri(const std::string& uri)
{
    return FromUris(std::vector<std::string>(1, uri));
}

std::unique_ptr<DataSetBase> DataSetIO::FromUris(const std::vector<std::string>& uris)
{
    if (uris.empty()) throw std::runtime_error{"[pbbam] dataset I/O ERROR: empty input URI list"};

    // create dataset(s) from URI(s)
    std::vector<std::unique_ptr<DataSetBase> > datasets;
    datasets.reserve(uris.size());
    for (const auto& uri : uris)
        datasets.emplace_back(DataSetFromUri(uri));
    assert(!datasets.empty());

    // if only 1, just return
    if (datasets.size() == 1) return std::unique_ptr<DataSetBase>(datasets.front().release());

    // else merge
    else {
        auto& result = datasets.at(0);
        for (size_t i = 1; i < datasets.size(); ++i) {
            const auto& next = datasets.at(i);
            *result += *next;
        }
        return std::move(result);
    }
}

std::unique_ptr<DataSetBase> DataSetIO::FromXmlString(const std::string& xml)
{
    if (xml.empty())
        throw std::runtime_error{"[pbbam] dataset I/O ERROR: cannot load from empty XML string"};
    std::istringstream s{xml};
    return XmlReader::FromStream(s);
}

void DataSetIO::ToFile(const std::unique_ptr<DataSetBase>& dataset, const std::string& fn,
                       DataSetPathMode pathMode)
{
    DataSetIO::ToFile(*dataset, fn, pathMode);
}

void DataSetIO::ToStream(const std::unique_ptr<DataSetBase>& dataset, std::ostream& out,
                         DataSetPathMode pathMode)
{
    DataSetIO::ToStream(*dataset, out, pathMode);
}

void DataSetIO::ToFile(DataSetBase& dataset, const std::string& fn, DataSetPathMode pathMode)
{
    std::ofstream out(fn);
    if (!out) throw DataSetFileException{fn, "could not open XML file for writing"};
    XmlWriter::ToStream(dataset, out, pathMode);
}

void DataSetIO::ToStream(DataSetBase& dataset, std::ostream& out, DataSetPathMode pathMode)
{
    XmlWriter::ToStream(dataset, out, pathMode);
}

}  // namespace BAM
}  // namespace PacBio
