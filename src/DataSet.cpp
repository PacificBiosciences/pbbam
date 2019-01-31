// File Description
/// \file DataSet.cpp
/// \brief Implements the DataSet class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/DataSet.h"

#include <algorithm>
#include <map>
#include <unordered_map>

#include <boost/algorithm/string.hpp>
#include <boost/icl/interval_set.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/optional.hpp>

#include "DataSetIO.h"
#include "FileUtils.h"
#include "TimeUtils.h"

#include "pbbam/DataSetTypes.h"
#include "pbbam/MakeUnique.h"
#include "pbbam/internal/DataSetBaseTypes.h"

namespace PacBio {
namespace BAM {
namespace {

const std::string defaultVersion{"4.0.0"};

void GetAllFiles(const ExternalResources& resources, std::vector<std::string>* result)
{
    for (const auto& resource : resources) {

        // store this resource's path
        result->push_back(resource.ResourceId());

        // store any child indices
        for (const auto& idx : resource.FileIndices())
            result->push_back(idx.ResourceId());

        // recurse into any other child resources
        GetAllFiles(resource.ExternalResources(), result);
    }
}

void InitDefaults(DataSet& ds)
{
    // provide default 'CreatedAt' & 'Version' attributes if not already present in XML

    if (ds.CreatedAt().empty()) ds.CreatedAt(TimeUtils::ToIso8601(TimeUtils::CurrentTime()));
    if (ds.Version().empty()) ds.Version(defaultVersion);
}

}  // anonymous

using internal::DataSetElement;

DataSet::DataSet() : DataSet(DataSet::GENERIC) { InitDefaults(*this); }

DataSet::DataSet(const DataSet::TypeEnum type) : path_(FileUtils::CurrentWorkingDirectory())
{
    switch (type) {
        case DataSet::GENERIC:
            d_ = std::make_unique<DataSetBase>();
            break;
        case DataSet::ALIGNMENT:
            d_ = std::make_unique<AlignmentSet>();
            break;
        case DataSet::BARCODE:
            d_ = std::make_unique<BarcodeSet>();
            break;
        case DataSet::CONSENSUS_ALIGNMENT:
            d_ = std::make_unique<ConsensusAlignmentSet>();
            break;
        case DataSet::CONSENSUS_READ:
            d_ = std::make_unique<ConsensusReadSet>();
            break;
        case DataSet::CONTIG:
            d_ = std::make_unique<ContigSet>();
            break;
        case DataSet::HDF_SUBREAD:
            d_ = std::make_unique<HdfSubreadSet>();
            break;
        case DataSet::REFERENCE:
            d_ = std::make_unique<ReferenceSet>();
            break;
        case DataSet::SUBREAD:
            d_ = std::make_unique<SubreadSet>();
            break;
        case DataSet::TRANSCRIPT:
            d_ = std::make_unique<TranscriptSet>();
            break;
        case DataSet::TRANSCRIPT_ALIGNMENT:
            d_ = std::make_unique<TranscriptAlignmentSet>();
            break;
        default:
            throw std::runtime_error{"unsupported dataset type"};
    }

    InitDefaults(*this);
}

DataSet::DataSet(const BamFile& bamFile)
    : d_(DataSetIO::FromUri(bamFile.Filename())), path_(FileUtils::CurrentWorkingDirectory())
{
    InitDefaults(*this);
}

DataSet::DataSet(const std::string& filename)
    : d_(DataSetIO::FromUri(filename)), path_(FileUtils::DirectoryName(filename))
{
    // for FOFN contents and raw BAM filenames, we can just use the current
    // directory as the starting path.
    //
    // (any relative paths in the FOFN have already been resolved)
    //
    if (boost::algorithm::iends_with(filename, ".fofn") ||
        boost::algorithm::iends_with(filename, ".bam") ||
        boost::algorithm::iends_with(filename, ".fasta") ||
        boost::algorithm::iends_with(filename, ".fa")) {
        path_ = FileUtils::CurrentWorkingDirectory();
    }
    InitDefaults(*this);
}

DataSet::DataSet(const std::vector<std::string>& filenames)
    : d_(DataSetIO::FromUris(filenames)), path_(FileUtils::CurrentWorkingDirectory())
{
    InitDefaults(*this);
}

DataSet::DataSet(const DataSet& other) : path_(other.path_)
{
    std::ostringstream out;
    DataSetIO::ToStream(other.d_, out);
    const std::string xml = out.str();
    d_ = DataSetIO::FromXmlString(xml);
}

DataSet& DataSet::operator=(const DataSet& other)
{
    if (this != &other) {
        std::ostringstream out;
        DataSetIO::ToStream(other.d_, out);
        const std::string xml = out.str();
        d_ = DataSetIO::FromXmlString(xml);
        path_ = other.path_;
    }
    return *this;
}

DataSet& DataSet::operator+=(const DataSet& other)
{
    *d_.get() += *other.d_.get();
    return *this;
}

std::vector<std::string> DataSet::AllFiles() const
{
    // get all files
    std::vector<std::string> result;
    GetAllFiles(ExternalResources(), &result);

    // resolve relative paths
    std::transform(result.begin(), result.end(), result.begin(),
                   [this](const std::string& fn) { return this->ResolvePath(fn); });
    return result;
}

std::vector<BamFile> DataSet::BamFiles() const
{
    const PacBio::BAM::ExternalResources& resources = ExternalResources();

    std::vector<BamFile> result;
    result.reserve(resources.Size());
    for (const ExternalResource& ext : resources) {

        // only bother resolving file path if this is a BAM file
        boost::iterator_range<std::string::const_iterator> bamFound =
            boost::algorithm::ifind_first(ext.MetaType(), "bam");
        if (!bamFound.empty()) {
            const std::string fn = ResolvePath(ext.ResourceId());
            result.emplace_back(fn);
        }
    }
    return result;
}

std::vector<std::string> DataSet::FastaFiles() const
{
    const PacBio::BAM::ExternalResources& resources = ExternalResources();

    std::vector<std::string> result;
    result.reserve(resources.Size());
    for (const ExternalResource& ext : resources) {

        // only bother resolving file path if this is a BAM file
        boost::iterator_range<std::string::const_iterator> fastaFound =
            boost::algorithm::ifind_first(ext.MetaType(), "fasta");
        if (!fastaFound.empty()) {
            const std::string fn = ResolvePath(ext.ResourceId());
            result.push_back(fn);
        }
    }
    return result;
}

DataSet DataSet::FromXml(const std::string& xml)
{
    DataSet result;
    result.d_ = DataSetIO::FromXmlString(xml);
    InitDefaults(result);
    return result;
}

const NamespaceRegistry& DataSet::Namespaces() const { return d_->Namespaces(); }

NamespaceRegistry& DataSet::Namespaces() { return d_->Namespaces(); }

DataSet::TypeEnum DataSet::NameToType(const std::string& typeName)
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
        lookup["TranscriptSet"] = DataSet::TRANSCRIPT;
        lookup["TranscriptAlignmentSet"] = DataSet::TRANSCRIPT_ALIGNMENT;
    }
    return lookup.at(typeName);  // throws if unknown typename
}

std::vector<std::string> DataSet::ResolvedResourceIds() const
{
    const PacBio::BAM::ExternalResources& resources = ExternalResources();

    std::vector<std::string> result;
    result.reserve(resources.Size());
    for (const ExternalResource& ext : resources) {
        result.push_back(ResolvePath(ext.ResourceId()));
    }
    return result;
}

std::string DataSet::ResolvePath(const std::string& originalPath) const
{
    return FileUtils::ResolvedFilePath(originalPath, path_);
}

void DataSet::Save(const std::string& outputFilename) const
{
    DataSetIO::ToFile(d_, outputFilename);
}

void DataSet::SaveToStream(std::ostream& out) const { DataSetIO::ToStream(d_, out); }

std::set<std::string> DataSet::SequencingChemistries() const
{
    const std::vector<BamFile> bamFiles{BamFiles()};

    std::set<std::string> result;
    for (const BamFile& bf : bamFiles) {
        if (!bf.IsPacBioBAM()) throw std::runtime_error{"only PacBio BAMs are supported"};
        const std::vector<ReadGroupInfo> readGroups{bf.Header().ReadGroups()};
        for (const ReadGroupInfo& rg : readGroups)
            result.insert(rg.SequencingChemistry());
    }
    return result;
}

std::vector<GenomicInterval> DataSet::GenomicIntervals() const
{
    // need to gather the contig lengths
    std::map<std::string, int32_t> contigLengths;
    for (const BamFile& b : BamFiles()) {
        const BamHeader& header = b.Header();
        const int32_t numContigs = header.NumSequences();
        for (int32_t i = 0; i < numContigs; ++i) {
            const std::string refName = header.SequenceName(i);
            const int32_t refLength = boost::lexical_cast<int32_t>(header.SequenceLength(i));

            const auto it = contigLengths.find(refName);
            if (it == contigLengths.cend())
                contigLengths.emplace(refName, refLength);
            else if (it->second != refLength)
                throw std::runtime_error{refName + " occurs twice with different lengths ('" +
                                         std::to_string(it->second) + "' and '" +
                                         std::to_string(refLength) + "')"};
        }
    }

    // with the lengths of all contigs known, we can build
    // the minimal interval set induced by the filters
    using intT = boost::icl::interval_set<int32_t>;
    using intInterval = intT::interval_type;

    std::map<std::string, intT> contigIntervals;
    int32_t numFilters = 0;

    for (const auto& xmlFilter : Filters()) {
        ++numFilters;
        boost::optional<std::string> contigName;

        intT intersectedInterval{intInterval{0, std::numeric_limits<int32_t>::max()}};

        for (const auto& xmlProperty : xmlFilter.Properties()) {
            const std::string XmlName = xmlProperty.Name();
            const std::string XmlOperator = xmlProperty.Operator();
            const std::string XmlValue = xmlProperty.Value();

            if ("rname" == XmlName) {
                if ("=" == XmlOperator) {
                    contigName = XmlValue;

                    const auto it = contigLengths.find(XmlValue);
                    if (it == contigLengths.cend())
                        throw std::runtime_error{"Could not find contig '" + XmlValue +
                                                 "' in BAM files"};
                    else
                        intersectedInterval &= intInterval(0, it->second);
                } else
                    throw std::runtime_error{
                        '\'' + XmlOperator +
                        "' is an unrecognized property operator, only '=' is recognized"};
            } else if ("tstart" == XmlName) {
                if ((XmlOperator != "<") && (XmlOperator != "<="))
                    throw std::runtime_error{"tstart only supports '<' and '<=' operators"};

                const int32_t end = boost::lexical_cast<int32_t>(XmlValue) + ("<=" == XmlOperator);
                intersectedInterval &= intInterval(0, end);
            } else if ("tend" == XmlName) {
                if ((XmlOperator != ">") && (XmlOperator != ">="))
                    throw std::runtime_error{"tend only supports '>' and '>=' operators"};

                const int32_t start =
                    boost::lexical_cast<int32_t>(XmlValue) - (">=" == XmlOperator);
                intersectedInterval &= intInterval(start, std::numeric_limits<int32_t>::max());
            } else
                throw std::runtime_error{'\'' + XmlName +
                                         "' is an unrecognized filter property name"};
        }

        if (contigName)
            contigIntervals[contigName.value()] |= intersectedInterval;
        else
            throw std::runtime_error{"Current filter does not have a valid 'rname' attribute"};
    }

    // extract all GenomicIntervals
    std::vector<GenomicInterval> result;
    if (numFilters) {
        // have some filters, only return regions passing filters
        for (const auto& contigs : contigIntervals) {
            const std::string& contigName = contigs.first;
            for (const auto& i : contigs.second) {
                // don't append empty intervals to the result
                if (boost::icl::length(i)) result.emplace_back(contigName, i.lower(), i.upper());
            }
        }
    } else {
        // no filters, return complete list of intervals
        for (const auto& contigs : contigLengths)
            result.emplace_back(contigs.first, 0, contigs.second);
    }

    return result;
}

std::string DataSet::TypeToName(const DataSet::TypeEnum& type)
{
    switch (type) {
        case DataSet::GENERIC:
            return "DataSet";
        case DataSet::ALIGNMENT:
            return "AlignmentSet";
        case DataSet::BARCODE:
            return "BarcodeSet";
        case DataSet::CONSENSUS_ALIGNMENT:
            return "ConsensusAlignmentSet";
        case DataSet::CONSENSUS_READ:
            return "ConsensusReadSet";
        case DataSet::CONTIG:
            return "ContigSet";
        case DataSet::HDF_SUBREAD:
            return "HdfSubreadSet";
        case DataSet::REFERENCE:
            return "ReferenceSet";
        case DataSet::SUBREAD:
            return "SubreadSet";
        case DataSet::TRANSCRIPT:
            return "TranscriptSet";
        case DataSet::TRANSCRIPT_ALIGNMENT:
            return "TranscriptAlignmentSet";
        default:
            throw std::runtime_error{"unsupported dataset type"};
    }
}

// Exposed timestamp utils

std::string CurrentTimestamp() { return TimeUtils::ToDataSetFormat(TimeUtils::CurrentTime()); }

std::string ToDataSetFormat(const std::chrono::system_clock::time_point& tp)
{
    return TimeUtils::ToDataSetFormat(tp);
}

std::string ToDataSetFormat(const time_t& t)
{
    return TimeUtils::ToDataSetFormat(std::chrono::system_clock::from_time_t(t));
}

std::string ToIso8601(const std::chrono::system_clock::time_point& tp)
{
    return TimeUtils::ToIso8601(tp);
}

std::string ToIso8601(const time_t& t)
{
    return TimeUtils::ToIso8601(std::chrono::system_clock::from_time_t(t));
}

}  // namespace BAM
}  // namespace PacBio
