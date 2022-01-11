#include "PbbamInternalConfig.h"

#include <pbbam/DataSet.h>

#include <pbbam/internal/DataSetBaseTypes.h>
#include "DataSetIO.h"
#include "DataSetUtils.h"
#include "FileUtils.h"
#include "TimeUtils.h"

#include <boost/algorithm/string.hpp>
#include <boost/icl/interval_set.hpp>
#include <boost/lexical_cast.hpp>

#include <algorithm>
#include <map>
#include <optional>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

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
        for (const auto& idx : resource.FileIndices()) {
            result->push_back(idx.ResourceId());
        }

        // recurse into any other child resources
        GetAllFiles(resource.ExternalResources(), result);
    }
}

}  // namespace

using internal::DataSetElement;

DataSet::DataSet() : DataSet(DataSet::GENERIC) {}

DataSet::DataSet(const DataSet::TypeEnum type)
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
            throw std::runtime_error{"[pbbam] dataset ERROR: encountered unsupported type"};
    }

    d_->Path(FileUtils::CurrentWorkingDirectory());
}

DataSet::DataSet(const BamFile& bamFile) : d_(DataSetIO::FromUri(bamFile.Filename()))
{
    d_->Path(FileUtils::CurrentWorkingDirectory());
}

DataSet::DataSet(const std::string& filename) : d_(DataSetIO::FromUri(filename))
{
    // for FOFN contents and raw BAM filenames, we can just use the current
    // directory as the starting path.
    //
    // (any relative paths in the FOFN have already been resolved)
    //
    if (boost::algorithm::iends_with(filename, ".fofn") ||
        boost::algorithm::iends_with(filename, ".bam") ||
        boost::algorithm::iends_with(filename, ".fasta") ||
        boost::algorithm::iends_with(filename, ".fsa") ||
        boost::algorithm::iends_with(filename, ".fa")) {
        d_->Path(FileUtils::CurrentWorkingDirectory());
    }

    else {
        if (boost::algorithm::iends_with(filename, ".xml")) {
            d_->FromInputXml(true);
        }
        d_->Path(FileUtils::DirectoryName(filename));
    }
}

DataSet::DataSet(const std::vector<std::string>& filenames) : d_(DataSetIO::FromUris(filenames))
{
    d_->Path(FileUtils::CurrentWorkingDirectory());
}

DataSet::DataSet(const DataSet& other)
{
    const bool otherFromXml = other.d_->FromInputXml();
    std::ostringstream out;
    DataSetIO::ToStream(other.d_, out, DataSetPathMode::ALLOW_RELATIVE);
    const std::string xml = out.str();
    d_ = DataSetIO::FromXmlString(xml);
    d_->Path(other.d_->Path());
    d_->FromInputXml(otherFromXml);
}

DataSet& DataSet::operator=(const DataSet& other)
{
    if (this != &other) {
        *this = DataSet{other};
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

const std::string& DataSet::Attribute(const std::string& name) const { return d_->Attribute(name); }

std::string& DataSet::Attribute(const std::string& name) { return d_->Attribute(name); }

DataSet& DataSet::Attribute(const std::string& name, const std::string& value)
{
    d_->Attribute(name, value);
    return *this;
}

std::vector<BamFile> DataSet::BamFiles() const
{
    const auto fns = BamFilenames();
    std::vector<BamFile> result;
    result.reserve(fns.size());
    for (const auto& fn : fns) {
        result.emplace_back(fn);
    }
    return result;
}

std::vector<std::string> DataSet::BamFilenames() const
{
    std::vector<std::string> result;
    const BAM::ExternalResources& resources = ExternalResources();

    result.reserve(resources.Size());
    for (const ExternalResource& ext : resources) {

        // only bother resolving file path if this is a BAM file
        boost::iterator_range<std::string::const_iterator> bamFound =
            boost::algorithm::ifind_first(ext.MetaType(), "bam");
        if (!bamFound.empty()) {
            const auto fn = ResolvePath(ext.ResourceId());
            result.emplace_back(fn);
        }
    }
    return result;
}

const std::string& DataSet::CreatedAt() const { return d_->CreatedAt(); }

std::string& DataSet::CreatedAt() { return d_->CreatedAt(); }

DataSet& DataSet::CreatedAt(const std::string& createdAt)
{
    d_->CreatedAt(createdAt);
    return *this;
}

const BAM::Extensions& DataSet::Extensions() const { return d_->Extensions(); }

BAM::Extensions& DataSet::Extensions() { return d_->Extensions(); }

DataSet& DataSet::Extensions(const BAM::Extensions& extensions)
{
    d_->Extensions(extensions);
    return *this;
}

const BAM::ExternalResources& DataSet::ExternalResources() const { return d_->ExternalResources(); }

BAM::ExternalResources& DataSet::ExternalResources() { return d_->ExternalResources(); }

DataSet& DataSet::ExternalResources(const BAM::ExternalResources& resources)
{
    d_->ExternalResources(resources);
    return *this;
}

std::vector<std::string> DataSet::FastaFiles() const
{
    const BAM::ExternalResources& resources = ExternalResources();

    std::vector<std::string> result;
    result.reserve(resources.Size());
    for (const ExternalResource& ext : resources) {

        // only bother resolving file path if this is a BAM file
        boost::iterator_range<std::string::const_iterator> fastaFound =
            boost::algorithm::ifind_first(ext.MetaType(), "fasta");
        if (!fastaFound.empty()) {
            const auto fn = ResolvePath(ext.ResourceId());
            result.push_back(fn);
        }
    }
    return result;
}

const BAM::Filters& DataSet::Filters() const { return d_->Filters(); }

BAM::Filters& DataSet::Filters() { return d_->Filters(); }

DataSet& DataSet::Filters(const BAM::Filters& filters)
{
    d_->Filters(filters);
    return *this;
}

const std::string& DataSet::Format() const { return d_->Format(); }

std::string& DataSet::Format() { return d_->Format(); }

DataSet& DataSet::Format(const std::string& format)
{
    d_->Format(format);
    return *this;
}

DataSet DataSet::FromXml(const std::string& xml)
{
    DataSet result;
    result.d_ = DataSetIO::FromXmlString(xml);
    result.d_->Path(FileUtils::DirectoryName(xml));
    result.d_->FromInputXml(true);
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
            if (it == contigLengths.cend()) {
                contigLengths.emplace(refName, refLength);
            } else if (it->second != refLength) {
                throw std::runtime_error{"[pbbam] dataset ERROR: " + refName +
                                         " occurs twice with different lengths ('" +
                                         std::to_string(it->second) + "' and '" +
                                         std::to_string(refLength) + "')"};
            }
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
        std::optional<std::string> contigName;

        intT intersectedInterval{intInterval{0, std::numeric_limits<int32_t>::max()}};

        for (const auto& xmlProperty : xmlFilter.Properties()) {
            const auto& XmlName = xmlProperty.Name();
            const auto& XmlOperator = xmlProperty.Operator();
            const auto& XmlValue = xmlProperty.Value();

            if ("rname" == XmlName) {
                if ("=" == XmlOperator) {
                    contigName = XmlValue;

                    const auto it = contigLengths.find(XmlValue);
                    if (it == contigLengths.cend()) {
                        throw std::runtime_error{"[pbbam] dataset ERROR: could not find contig '" +
                                                 XmlValue + "' in BAM files"};
                    } else {
                        intersectedInterval &= intInterval(0, it->second);
                    }
                } else {
                    throw std::runtime_error{
                        "[pbbam] dataset ERROR: '" + XmlOperator +
                        "' is an unrecognized property operator, only '=' is recognized"};
                }
            } else if ("tstart" == XmlName) {
                if ((XmlOperator != "<") && (XmlOperator != "<=")) {
                    throw std::runtime_error{
                        "[pbbam] dataset ERROR: 'tstart' only supports '<' and '<=' operators"};
                }

                const int32_t end = boost::lexical_cast<int32_t>(XmlValue) + ("<=" == XmlOperator);
                intersectedInterval &= intInterval(0, end);
            } else if ("tend" == XmlName) {
                if ((XmlOperator != ">") && (XmlOperator != ">=")) {
                    throw std::runtime_error{
                        "[pbbam] dataset ERROR: 'tend' only supports '>' and '>=' operators"};
                }

                const int32_t start =
                    boost::lexical_cast<int32_t>(XmlValue) - (">=" == XmlOperator);
                intersectedInterval &= intInterval(start, std::numeric_limits<int32_t>::max());
            } else {
                throw std::runtime_error{"[pbbam] dataset ERROR: '" + XmlName +
                                         "' is an unrecognized filter property name"};
            }
        }

        if (contigName) {
            contigIntervals[contigName.value()] |= intersectedInterval;
        } else {
            throw std::runtime_error{
                "[pbbam] dataset ERROR: current filter does not have a valid 'rname' attribute"};
        }
    }

    // extract all GenomicIntervals
    std::vector<GenomicInterval> result;
    if (numFilters) {
        // have some filters, only return regions passing filters
        for (const auto& contigs : contigIntervals) {
            const auto& contigName = contigs.first;
            for (const auto& i : contigs.second) {
                // don't append empty intervals to the result
                if (boost::icl::length(i)) {
                    result.emplace_back(contigName, i.lower(), i.upper());
                }
            }
        }
    } else {
        // no filters, return complete list of intervals
        for (const auto& contigs : contigLengths) {
            result.emplace_back(contigs.first, 0, contigs.second);
        }
    }

    return result;
}

BamHeader DataSet::MergedHeader() const { return BamHeader{*this}; }

const BAM::DataSetMetadata& DataSet::Metadata() const { return d_->Metadata(); }

BAM::DataSetMetadata& DataSet::Metadata() { return d_->Metadata(); }

DataSet& DataSet::Metadata(const BAM::DataSetMetadata& metadata)
{
    d_->Metadata(metadata);
    return *this;
}

const std::string& DataSet::MetaType() const { return d_->MetaType(); }

std::string& DataSet::MetaType() { return d_->MetaType(); }

DataSet& DataSet::MetaType(const std::string& metatype)
{
    d_->MetaType(metatype);
    return *this;
}

const std::string& DataSet::ModifiedAt() const { return d_->ModifiedAt(); }

std::string& DataSet::ModifiedAt() { return d_->ModifiedAt(); }

DataSet& DataSet::ModifiedAt(const std::string& modifiedAt)
{
    d_->ModifiedAt(modifiedAt);
    return *this;
}

const std::string& DataSet::Name() const { return d_->Name(); }

std::string& DataSet::Name() { return d_->Name(); }

DataSet& DataSet::Name(const std::string& name)
{
    d_->Name(name);
    return *this;
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

const std::string& DataSet::Path() const { return d_->Path(); }

std::vector<std::string> DataSet::ResolvedResourceIds() const
{
    const BAM::ExternalResources& resources = ExternalResources();

    std::vector<std::string> result;
    result.reserve(resources.Size());
    for (const ExternalResource& ext : resources) {
        result.push_back(ResolvePath(ext.ResourceId()));
    }
    return result;
}

std::string DataSet::ResolvePath(const std::string& originalPath) const
{
    return FileUtils::ResolvedFilePath(originalPath, d_->Path());
}

const std::string& DataSet::ResourceId() const { return d_->ResourceId(); }

std::string& DataSet::ResourceId() { return d_->ResourceId(); }

DataSet& DataSet::ResourceId(const std::string& resourceId)
{
    d_->ResourceId(resourceId);
    return *this;
}

void DataSet::Save(const std::string& outputFilename, DataSetPathMode pathMode) const
{
    DataSetIO::ToFile(d_, outputFilename, pathMode);
}

void DataSet::SaveToStream(std::ostream& out, DataSetPathMode pathMode) const
{
    DataSetIO::ToStream(d_, out, pathMode);
}

std::set<std::string> DataSet::SequencingChemistries() const
{
    std::set<std::string> result;
    for (const auto& bf : BamFiles()) {
        if (!bf.IsPacBioBAM()) {
            throw std::runtime_error{
                "[pbbam] dataset ERROR: only PacBio BAMs are supported for fetching chemistry "
                "info"};
        }
        const std::vector<ReadGroupInfo> readGroups{bf.Header().ReadGroups()};
        for (const auto& rg : readGroups) {
            result.insert(rg.SequencingChemistry());
        }
    }
    return result;
}

std::set<std::string> DataSet::Samples() const
{
    std::set<std::string> result;
    for (const auto& bf : BamFiles()) {
        const std::vector<ReadGroupInfo> readGroups{bf.Header().ReadGroups()};
        for (const auto& rg : readGroups) {
            result.insert(rg.Sample());
        }
    }
    return result;
}

const BAM::SubDataSets& DataSet::SubDataSets() const { return d_->SubDataSets(); }

BAM::SubDataSets& DataSet::SubDataSets() { return d_->SubDataSets(); }

DataSet& DataSet::SubDataSets(const BAM::SubDataSets& subdatasets)
{
    d_->SubDataSets(subdatasets);
    return *this;
}

const BAM::SupplementalResources& DataSet::SupplementalResources() const
{
    return d_->SupplementalResources();
}

BAM::SupplementalResources& DataSet::SupplementalResources() { return d_->SupplementalResources(); }

DataSet& DataSet::SupplementalResources(const BAM::SupplementalResources& resources)
{
    d_->SupplementalResources(resources);
    return *this;
}

const std::string& DataSet::Tags() const { return d_->Tags(); }

std::string& DataSet::Tags() { return d_->Tags(); }

DataSet& DataSet::Tags(const std::string& tags)
{
    d_->Tags(tags);
    return *this;
}

const std::string& DataSet::TimeStampedName() const { return d_->TimeStampedName(); }

std::string& DataSet::TimeStampedName() { return d_->TimeStampedName(); }

DataSet& DataSet::TimeStampedName(const std::string& timeStampedName)
{
    d_->TimeStampedName(timeStampedName);
    return *this;
}

PacBio::BAM::DataSet::TypeEnum DataSet::Type() const { return DataSet::NameToType(TypeName()); }

DataSet& DataSet::Type(const DataSet::TypeEnum type)
{
    d_->Label(DataSet::TypeToName(type));
    return *this;
}

std::string DataSet::TypeName() const { return d_->LocalNameLabel().to_string(); }

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
            throw std::runtime_error{"[pbbam] dataset ERROR: encountered unsupported dataset type"};
    }
}

const std::string& DataSet::UniqueId() const { return d_->UniqueId(); }

std::string& DataSet::UniqueId() { return d_->UniqueId(); }

DataSet& DataSet::UniqueId(const std::string& uuid)
{
    d_->UniqueId(uuid);
    return *this;
}

const std::string& DataSet::Version() const { return d_->Version(); }

std::string& DataSet::Version() { return d_->Version(); }

DataSet& DataSet::Version(const std::string& version)
{
    d_->Version(version);
    return *this;
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
