// File Description
/// \file DataSetTypes.cpp
/// \brief Implementations for the public DataSet component classes.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/DataSetTypes.h"

#include <cstddef>
#include <set>
#include <unordered_map>

#include "pbbam/Unused.h"
#include "pbbam/internal/DataSetBaseTypes.h"

#include "DataSetIO.h"
#include "DataSetUtils.h"
#include "FileUtils.h"
#include "TimeUtils.h"

namespace {

// clang-format off
using ElementType = PacBio::BAM::XmlElementType;
const std::unordered_map<std::string, PacBio::BAM::XmlElementType> elementTypeLookup
{
    {"DataSetMetadata",        ElementType::DATASET_METADATA},
    {"ExtensionElement",       ElementType::EXTENSION},
    {"Extensions",             ElementType::EXTENSIONS},
    {"ExternalResource",       ElementType::EXTERNAL_RESOURCE},
    {"ExternalResources",      ElementType::EXTERNAL_RESOURCES},
    {"FileIndex",              ElementType::FILE_INDEX},
    {"FileIndices",            ElementType::FILE_INDICES},
    {"Filter",                 ElementType::FILTER},
    {"Filters",                ElementType::FILTERS},
    {"ParentTool",             ElementType::PARENT_TOOL},
    {"Property",               ElementType::PROPERTY},
    {"Properties",             ElementType::PROPERTIES},
    {"Provenance",             ElementType::PROVENANCE},
    {"AlignmentSet",           ElementType::ALIGNMENT_SET},
    {"BarcodeSet",             ElementType::BARCODE_SET},
    {"ConsensusAlignmentSet",  ElementType::CONSENSUS_ALIGNMENT_SET},
    {"ConsensusReadSet",       ElementType::CONSENSUS_READ_SET},
    {"ContigSet",              ElementType::CONTIG_SET},
    {"HdfSubreadSet",          ElementType::HDF_SUBREAD_SET},
    {"ReferenceSet",           ElementType::REFERENCE_SET},
    {"SubreadSet",             ElementType::SUBREAD_SET},
    {"TranscriptSet",          ElementType::TRANSCRIPT_SET},
    {"TranscriptAlignmentSet", ElementType::TRANSCRIPT_ALIGNMENT_SET},
    {"DataSets",               ElementType::SUBDATASETS},
    {"DataSet",                ElementType::GENERIC_DATASET}
};
// clang-format on

}  // anonymous

namespace PacBio {
namespace BAM {

// -------------------
// AlignmentSet
// -------------------

AlignmentSet::AlignmentSet()
    : DataSetBase("PacBio.DataSet.AlignmentSet", "AlignmentSet", XsdType::DATASETS)
{
}

AlignmentSet::AlignmentSet(const internal::FromInputXml& fromInputXml)
    : DataSetBase("", "AlignmentSet", fromInputXml, XsdType::DATASETS)
{
}

// -------------------
// BarcodeSet
// -------------------

BarcodeSet::BarcodeSet() : DataSetBase("PacBio.DataSet.BarcodeSet", "BarcodeSet", XsdType::DATASETS)
{
}

BarcodeSet::BarcodeSet(const internal::FromInputXml& fromInputXml)
    : DataSetBase("", "BarcodeSet", fromInputXml, XsdType::DATASETS)
{
}

// -----------------------
// ConsensusAlignmentSet
// -----------------------

ConsensusAlignmentSet::ConsensusAlignmentSet()
    : DataSetBase("PacBio.DataSet.ConsensusAlignmentSet", "ConsensusAlignmentSet",
                  XsdType::DATASETS)
{
}

ConsensusAlignmentSet::ConsensusAlignmentSet(const internal::FromInputXml& fromInputXml)
    : DataSetBase("", "ConsensusAlignmentSet", fromInputXml, XsdType::DATASETS)
{
}

// -------------------
// ConsensusReadSet
// -------------------

ConsensusReadSet::ConsensusReadSet()
    : DataSetBase("PacBio.DataSet.ConsensusReadSet", "ConsensusReadSet", XsdType::DATASETS)
{
}

ConsensusReadSet::ConsensusReadSet(const internal::FromInputXml& fromInputXml)
    : DataSetBase("", "ConsensusReadSet", fromInputXml, XsdType::DATASETS)
{
}

// -------------------
// ContigSet
// -------------------

ContigSet::ContigSet() : DataSetBase("PacBio.DataSet.ContigSet", "ContigSet", XsdType::DATASETS) {}

ContigSet::ContigSet(const internal::FromInputXml& fromInputXml)
    : DataSetBase("", "ContigSet", fromInputXml, XsdType::DATASETS)
{
}

// -------------------
// DataSetBase
// -------------------

DataSetBase::DataSetBase()
    : StrictEntityType("PacBio.DataSet.DataSet", "DataSet", XsdType::DATASETS)
{
}

DataSetBase::DataSetBase(const internal::FromInputXml& fromInputXml)
    : StrictEntityType("", "DataSet", fromInputXml, XsdType::DATASETS)
{
}

DataSetBase::DataSetBase(const std::string& metatype, const std::string& label, const XsdType& xsd)
    : StrictEntityType(metatype, label, xsd)
{
}

DataSetBase::DataSetBase(const std::string& metatype, const std::string& label,
                         const internal::FromInputXml& fromInputXml, const XsdType& xsd)
    : StrictEntityType(metatype, label, fromInputXml, xsd)
{
}

const PacBio::BAM::ExternalResources& DataSetBase::ExternalResources() const
{
    return Child<PacBio::BAM::ExternalResources>("ExternalResources");
}

PacBio::BAM::ExternalResources& DataSetBase::ExternalResources()
{
    if (!HasChild("ExternalResources")) AddChild(PacBio::BAM::ExternalResources());
    auto& c = Child<PacBio::BAM::ExternalResources>("ExternalResources");
    return c;
}

DataSetBase& DataSetBase::ExternalResources(const PacBio::BAM::ExternalResources& resources)
{
    ExternalResources() = resources;
    return *this;
}

DEFINE_ACCESSORS(DataSetBase, Filters, Filters)

DataSetBase& DataSetBase::Filters(const PacBio::BAM::Filters& filters)
{
    Filters() = filters;
    return *this;
}

DEFINE_ACCESSORS(DataSetBase, DataSetMetadata, Metadata)

DataSetBase& DataSetBase::Metadata(const PacBio::BAM::DataSetMetadata& metadata)
{
    Metadata() = metadata;
    return *this;
}

const PacBio::BAM::SubDataSets& DataSetBase::SubDataSets() const
{
    try {
        return Child<PacBio::BAM::SubDataSets>("DataSets");
    } catch (std::exception&) {
        return internal::NullObject<PacBio::BAM::SubDataSets>();
    }
}

PacBio::BAM::SubDataSets& DataSetBase::SubDataSets()
{
    if (!HasChild("DataSets")) AddChild(internal::NullObject<PacBio::BAM::SubDataSets>());
    return Child<PacBio::BAM::SubDataSets>("DataSets");
}

DataSetBase& DataSetBase::SubDataSets(const PacBio::BAM::SubDataSets& subdatasets)
{
    SubDataSets() = subdatasets;
    return *this;
}

DataSetBase* DataSetBase::DeepCopy() const
{
    auto* copyDataset = new DataSetElement(*this);
    auto* result = static_cast<DataSetBase*>(copyDataset);
    result->registry_ = registry_;
    return result;
}

DataSetBase& DataSetBase::operator+=(const DataSetBase& other)
{
    // must be same dataset types (or 'other' must be generic)
    if (other.LocalNameLabel() != LocalNameLabel() && other.LocalNameLabel() != "DataSet")
        throw std::runtime_error{"cannot merge different dataset types"};

    // check object metadata
    Metadata() += other.Metadata();
    ExternalResources() += other.ExternalResources();
    Filters() += other.Filters();
    SubDataSets() += other;

    return *this;
}

std::shared_ptr<DataSetBase> DataSetBase::Create(const std::string& typeName)
{
    if (typeName == std::string("DataSet")) return std::make_shared<DataSetBase>();
    if (typeName == std::string("SubreadSet")) return std::make_shared<SubreadSet>();
    if (typeName == std::string("AlignmentSet")) return std::make_shared<AlignmentSet>();
    if (typeName == std::string("BarcodeSet")) return std::make_shared<BarcodeSet>();
    if (typeName == std::string("ConsensusAlignmentSet"))
        return std::make_shared<ConsensusAlignmentSet>();
    if (typeName == std::string("ConsensusReadSet")) return std::make_shared<ConsensusReadSet>();
    if (typeName == std::string("ContigSet")) return std::make_shared<ContigSet>();
    if (typeName == std::string("HdfSubreadSet")) return std::make_shared<HdfSubreadSet>();
    if (typeName == std::string("ReferenceSet")) return std::make_shared<ReferenceSet>();
    if (typeName == std::string("TranscriptSet")) return std::make_shared<TranscriptSet>();
    if (typeName == std::string("TranscriptAlignmentSet"))
        return std::make_shared<TranscriptAlignmentSet>();

    // unknown typename
    throw std::runtime_error{"unsupported dataset type: " + typeName};
}

std::shared_ptr<DataSetBase> DataSetBase::Create(const std::string& typeName,
                                                 const internal::FromInputXml& fromInputXml)
{
    if (typeName == std::string("DataSet")) return std::make_shared<DataSetBase>(fromInputXml);
    if (typeName == std::string("SubreadSet")) return std::make_shared<SubreadSet>(fromInputXml);
    if (typeName == std::string("AlignmentSet"))
        return std::make_shared<AlignmentSet>(fromInputXml);
    if (typeName == std::string("BarcodeSet")) return std::make_shared<BarcodeSet>(fromInputXml);
    if (typeName == std::string("ConsensusAlignmentSet"))
        return std::make_shared<ConsensusAlignmentSet>(fromInputXml);
    if (typeName == std::string("ConsensusReadSet"))
        return std::make_shared<ConsensusReadSet>(fromInputXml);
    if (typeName == std::string("ContigSet")) return std::make_shared<ContigSet>(fromInputXml);
    if (typeName == std::string("HdfSubreadSet"))
        return std::make_shared<HdfSubreadSet>(fromInputXml);
    if (typeName == std::string("ReferenceSet"))
        return std::make_shared<ReferenceSet>(fromInputXml);
    if (typeName == std::string("TranscriptSet"))
        return std::make_shared<TranscriptSet>(fromInputXml);
    if (typeName == std::string("TranscriptAlignmentSet"))
        return std::make_shared<TranscriptAlignmentSet>(fromInputXml);

    // unknown typename
    throw std::runtime_error{"unsupported dataset type: " + typeName};
}

void DataSetBase::Save(const std::string& outputFilename)
{
    DataSetIO::ToFile(*this, outputFilename);
}

void DataSetBase::SaveToStream(std::ostream& out) { DataSetIO::ToStream(*this, out); }

// -------------------
// DataSetMetadata
// -------------------

DataSetMetadata::DataSetMetadata() : DataSetElement("DataSetMetadata", XsdType::DATASETS) {}

DataSetMetadata::DataSetMetadata(const internal::FromInputXml& fromInputXml)
    : DataSetElement("", fromInputXml, XsdType::DATASETS)
{
}

DataSetMetadata::DataSetMetadata(const std::string& numRecords, const std::string& totalLength)
    : DataSetElement("DataSetMetadata", XsdType::DATASETS)
{
    TotalLength(totalLength);
    NumRecords(numRecords);
}

DataSetMetadata::DataSetMetadata(const std::string& numRecords, const std::string& totalLength,
                                 const internal::FromInputXml& fromInputXml)
    : DataSetElement("", fromInputXml, XsdType::DATASETS)
{
    TotalLength(totalLength);
    NumRecords(numRecords);
}

DEFINE_ACCESSORS(DataSetMetadata, Provenance, Provenance)

DataSetMetadata& DataSetMetadata::Provenance(const PacBio::BAM::Provenance& provenance)
{
    Provenance() = provenance;
    return *this;
}

DataSetMetadata& DataSetMetadata::operator+=(const DataSetMetadata& other)
{
    TotalLength() = TotalLength() + other.TotalLength();
    NumRecords() = NumRecords() + other.NumRecords();
    // merge add'l
    return *this;
}

// -------------------
// ExtensionElement
// -------------------

ExtensionElement::ExtensionElement() : DataSetElement("ExtensionElement", XsdType::BASE_DATA_MODEL)
{
}

ExtensionElement::ExtensionElement(const internal::FromInputXml& fromInputXml)
    : DataSetElement("", fromInputXml, XsdType::BASE_DATA_MODEL)
{
}

// -------------------
// Extensions
// -------------------

Extensions::Extensions() : DataSetElement("Extensions", XsdType::BASE_DATA_MODEL) {}

Extensions::Extensions(const internal::FromInputXml& fromInputXml)
    : DataSetElement("", fromInputXml, XsdType::BASE_DATA_MODEL)
{
}

// -------------------
// ExternalResource
// -------------------

ExternalResource::ExternalResource(const BamFile& bamFile)
    : IndexedDataType("PacBio.SubreadFile.SubreadBamFile", bamFile.Filename(), "ExternalResource",
                      XsdType::BASE_DATA_MODEL)
{
}

ExternalResource::ExternalResource(const std::string& metatype, const std::string& filename)
    : IndexedDataType(metatype, filename, "ExternalResource", XsdType::BASE_DATA_MODEL)
{
}

ExternalResource::ExternalResource(const std::string& metatype, const std::string& filename,
                                   const internal::FromInputXml& fromInputXml)
    : IndexedDataType("", filename, "ExternalResource", fromInputXml, XsdType::BASE_DATA_MODEL)
{
    UNUSED(metatype);
}

DEFINE_ACCESSORS(ExternalResource, ExternalResources, ExternalResources)

ExternalResource& ExternalResource::ExternalResources(
    const PacBio::BAM::ExternalResources& resources)
{
    ExternalResources() = resources;
    return *this;
}

BamFile ExternalResource::ToBamFile() const { return BamFile(ResourceId()); }

// -------------------
// ExternalResources
// -------------------

ExternalResources::ExternalResources()
    : DataSetElement("ExternalResources", XsdType::BASE_DATA_MODEL)
{
}

ExternalResources::ExternalResources(const internal::FromInputXml& fromInputXml)
    : DataSetElement("", fromInputXml, XsdType::BASE_DATA_MODEL)
{
}

ExternalResources& ExternalResources::operator+=(const ExternalResources& other)
{
    // only keep unique resource ids
    std::set<std::string> myResourceIds;
    for (size_t i = 0; i < NumChildren(); ++i) {
        const ExternalResource& resource = this->operator[](i);
        myResourceIds.insert(resource.ResourceId());
    }

    std::vector<size_t> newResourceIndices;
    const size_t numOtherResourceIds = other.Size();
    for (size_t i = 0; i < numOtherResourceIds; ++i) {
        const std::string& resourceId = other[i].ResourceId();
        auto found = myResourceIds.find(resourceId);
        if (found == myResourceIds.cend()) newResourceIndices.push_back(i);
    }

    for (size_t index : newResourceIndices)
        Add(other[index]);
    return *this;
}

void ExternalResources::Add(const ExternalResource& ext)
{
    // disallow external resources w/ duplicate ResourceIds
    std::set<std::string> myResourceIds;
    for (size_t i = 0; i < NumChildren(); ++i) {
        const ExternalResource& resource = this->operator[](i);
        myResourceIds.insert(resource.ResourceId());
    }

    if (myResourceIds.find(ext.ResourceId()) == myResourceIds.cend()) AddChild(ext);
}

std::vector<BamFile> ExternalResources::BamFiles() const
{

    std::vector<BamFile> result;
    const int numResources = Size();
    result.reserve(numResources);
    for (const ExternalResource& ext : *this)
        result.push_back(ext.ToBamFile());
    return result;
}

void ExternalResources::Remove(const ExternalResource& ext) { RemoveChild(ext); }

// -------------------
// FileIndex
// -------------------

FileIndex::FileIndex(const std::string& metatype, const std::string& filename)
    : InputOutputDataType(metatype, filename, "FileIndex", XsdType::BASE_DATA_MODEL)
{
}

FileIndex::FileIndex(const std::string& metatype, const std::string& filename,
                     const internal::FromInputXml& fromInputXml)
    : InputOutputDataType("", filename, "FileIndex", fromInputXml, XsdType::BASE_DATA_MODEL)
{
    UNUSED(metatype);
}

// -------------------
// FileIndices
// -------------------

FileIndices::FileIndices() : DataSetElement("FileIndices", XsdType::BASE_DATA_MODEL) {}

FileIndices::FileIndices(const internal::FromInputXml& fromInputXml)
    : DataSetElement("", fromInputXml, XsdType::BASE_DATA_MODEL)
{
}

void FileIndices::Add(const FileIndex& index) { AddChild(index); }

void FileIndices::Remove(const FileIndex& index) { RemoveChild(index); }

// -------------------
// Filter
// -------------------

Filter::Filter() : DataSetElement("Filter", XsdType::DATASETS) {}

Filter::Filter(const internal::FromInputXml& fromInputXml)
    : DataSetElement("Filter", fromInputXml, XsdType::DATASETS)
{
}

DEFINE_ACCESSORS(Filter, Properties, Properties)

Filter& Filter::Properties(const PacBio::BAM::Properties& properties)
{
    Properties() = properties;
    return *this;
}

// -------------------
// Filters
// -------------------

Filters::Filters() : DataSetElement("Filters", XsdType::DATASETS) {}

Filters::Filters(const internal::FromInputXml& fromInputXml)
    : DataSetElement("", fromInputXml, XsdType::DATASETS)
{
}

Filters& Filters::operator+=(const Filters& other)
{
    for (auto& newFilter : other)
        AddChild(newFilter);
    return *this;
}

void Filters::Add(const Filter& filter) { AddChild(filter); }

void Filters::Remove(const Filter& filter) { RemoveChild(filter); }

// -------------------
// HdfSubreadSet
// -------------------

HdfSubreadSet::HdfSubreadSet()
    : DataSetBase("PacBio.DataSet.HdfSubreadSet", "HdfSubreadSet", XsdType::DATASETS)
{
}

HdfSubreadSet::HdfSubreadSet(const internal::FromInputXml& fromInputXml)
    : DataSetBase("", "HdfSubreadSet", fromInputXml, XsdType::DATASETS)
{
}

// -------------------
// ParentTool
// -------------------

ParentTool::ParentTool() : BaseEntityType("ParentTool", XsdType::DATASETS) {}

ParentTool::ParentTool(const internal::FromInputXml& fromInputXml)
    : BaseEntityType("", fromInputXml, XsdType::DATASETS)
{
}

// -------------------
// Properties
// -------------------

Properties::Properties() : DataSetElement("Properties", XsdType::BASE_DATA_MODEL) {}

Properties::Properties(const internal::FromInputXml& fromInputXml)
    : DataSetElement("", fromInputXml, XsdType::BASE_DATA_MODEL)
{
}

void Properties::Add(const Property& property) { AddChild(property); }

void Properties::Remove(const Property& property) { RemoveChild(property); }

// -------------------
// Property
// -------------------

Property::Property(const std::string& name, const std::string& value, const std::string& op)
    : DataSetElement("Property", XsdType::BASE_DATA_MODEL)
{
    Name(name);
    Value(value);
    Operator(op);
}

Property::Property(const std::string& name, const std::string& value, const std::string& op,
                   const internal::FromInputXml& fromInputXml)
    : DataSetElement("", fromInputXml, XsdType::BASE_DATA_MODEL)
{
    Name(name);
    Value(value);
    Operator(op);
}

// -------------------
// Provenance
// -------------------

Provenance::Provenance() : DataSetElement("Provenance", XsdType::DATASETS) {}

Provenance::Provenance(const internal::FromInputXml& fromInputXml)
    : DataSetElement("", fromInputXml, XsdType::DATASETS)
{
}

DEFINE_ACCESSORS(Provenance, ParentTool, ParentTool)

// -------------------
// ReferenceSet
// -------------------

ReferenceSet::ReferenceSet()
    : DataSetBase("PacBio.DataSet.ReferenceSet", "ReferenceSet", XsdType::DATASETS)
{
}

ReferenceSet::ReferenceSet(const internal::FromInputXml& fromInputXml)
    : DataSetBase("", "ReferenceSet", fromInputXml, XsdType::DATASETS)
{
}

// -------------------
// SubDataSets
// -------------------

SubDataSets::SubDataSets() : DataSetElement("DataSets", XsdType::DATASETS) {}

SubDataSets::SubDataSets(const internal::FromInputXml& fromInputXml)
    : DataSetElement("", fromInputXml, XsdType::DATASETS)
{
}

SubDataSets& SubDataSets::operator+=(const DataSetBase& other)
{
    AddChild(other);
    return *this;
}

SubDataSets& SubDataSets::operator+=(const SubDataSets& other)
{
    for (auto& newSubDataset : other)
        AddChild(newSubDataset);
    return *this;
}

void SubDataSets::Add(const DataSetBase& subdataset) { AddChild(subdataset); }

void SubDataSets::Remove(const DataSetBase& subdataset) { RemoveChild(subdataset); }

// -------------------
// SubreadSet
// -------------------

SubreadSet::SubreadSet() : DataSetBase("PacBio.DataSet.SubreadSet", "SubreadSet", XsdType::DATASETS)
{
}

SubreadSet::SubreadSet(const internal::FromInputXml& fromInputXml)
    : DataSetBase("", "SubreadSet", fromInputXml, XsdType::DATASETS)
{
}

// -------------------
// TranscriptSet
// -------------------

TranscriptSet::TranscriptSet()
    : DataSetBase("PacBio.DataSet.TranscriptSet", "TranscriptSet", XsdType::DATASETS)
{
}

TranscriptSet::TranscriptSet(const internal::FromInputXml& fromInputXml)
    : DataSetBase("", "TranscriptSet", fromInputXml, XsdType::DATASETS)
{
}

// -------------------
// TranscriptAlignmentSet
// -------------------

TranscriptAlignmentSet::TranscriptAlignmentSet()
    : DataSetBase("PacBio.DataSet.TranscriptAlignmentSet", "TranscriptAlignmentSet",
                  XsdType::DATASETS)
{
}

TranscriptAlignmentSet::TranscriptAlignmentSet(const internal::FromInputXml& fromInputXml)
    : DataSetBase("", "TranscriptAlignmentSet", fromInputXml, XsdType::DATASETS)
{
}

XmlElementType ElementTypeFromName(const std::string& name)
{
    const auto found = elementTypeLookup.find(name);
    if (found == elementTypeLookup.cend()) return XmlElementType::GENERIC_ELEMENT;
    return found->second;
}

}  // namespace BAM
}  // namespace PacBio
