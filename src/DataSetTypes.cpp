// File Description
/// \file DataSetTypes.cpp
/// \brief Implementations for the public DataSet component classes.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/DataSetTypes.h"

#include <cstddef>
#include <set>

#include "DataSetUtils.h"
#include "FileUtils.h"
#include "TimeUtils.h"
#include "pbbam/internal/DataSetBaseTypes.h"

namespace PacBio {
namespace BAM {

// -------------------
// AlignmentSet
// -------------------

AlignmentSet::AlignmentSet()
    : DataSetBase("PacBio.DataSet.AlignmentSet", "AlignmentSet", XsdType::DATASETS)
{
}

// -------------------
// BarcodeSet
// -------------------

BarcodeSet::BarcodeSet() : DataSetBase("PacBio.DataSet.BarcodeSet", "BarcodeSet", XsdType::DATASETS)
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

// -------------------
// ConsensusReadSet
// -------------------

ConsensusReadSet::ConsensusReadSet()
    : DataSetBase("PacBio.DataSet.ConsensusReadSet", "ConsensusReadSet", XsdType::DATASETS)
{
}

// -------------------
// ContigSet
// -------------------

ContigSet::ContigSet() : DataSetBase("PacBio.DataSet.ContigSet", "ContigSet", XsdType::DATASETS) {}

// -------------------
// DataSetBase
// -------------------

DataSetBase::DataSetBase()
    : StrictEntityType("PacBio.DataSet.DataSet", "DataSet", XsdType::DATASETS)
{
}

DataSetBase::DataSetBase(const std::string& metatype, const std::string& label, const XsdType& xsd)
    : StrictEntityType(metatype, label, xsd)
{
}

DEFINE_ACCESSORS(DataSetBase, ExternalResources, ExternalResources)

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

    // check filter match
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

    // unknown typename
    throw std::runtime_error{"unsupported dataset type"};
}

// -------------------
// DataSetMetadata
// -------------------

DataSetMetadata::DataSetMetadata(const std::string& numRecords, const std::string& totalLength)
    : DataSetElement("DataSetMetadata", XsdType::DATASETS)
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

// -------------------
// Extensions
// -------------------

Extensions::Extensions()
    : DataSetListElement<ExtensionElement>("Extensions", XsdType::BASE_DATA_MODEL)
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
    : DataSetListElement<ExternalResource>("ExternalResources", XsdType::BASE_DATA_MODEL)
{
}

ExternalResources& ExternalResources::operator+=(const ExternalResources& other)
{
    // only keep unique resource ids

    std::set<std::string> myResourceIds;
    for (size_t i = 0; i < Size(); ++i) {
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
    for (size_t i = 0; i < Size(); ++i) {
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

// -------------------
// FileIndices
// -------------------

FileIndices::FileIndices() : DataSetListElement<FileIndex>("FileIndices", XsdType::BASE_DATA_MODEL)
{
}

void FileIndices::Add(const FileIndex& index) { AddChild(index); }

void FileIndices::Remove(const FileIndex& index) { RemoveChild(index); }

// -------------------
// Filter
// -------------------

Filter::Filter() : DataSetElement("Filter", XsdType::DATASETS) {}

DEFINE_ACCESSORS(Filter, Properties, Properties)

Filter& Filter::Properties(const PacBio::BAM::Properties& properties)
{
    Properties() = properties;
    return *this;
}

// -------------------
// Filters
// -------------------

Filters::Filters() : DataSetListElement<Filter>("Filters", XsdType::DATASETS) {}

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

// -------------------
// ParentTool
// -------------------

ParentTool::ParentTool() : BaseEntityType("ParentTool", XsdType::DATASETS) {}

// -------------------
// Properties
// -------------------

Properties::Properties() : DataSetListElement<Property>("Properties", XsdType::BASE_DATA_MODEL) {}

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

// -------------------
// Provenance
// -------------------

Provenance::Provenance() : DataSetElement("Provenance", XsdType::DATASETS) {}

DEFINE_ACCESSORS(Provenance, ParentTool, ParentTool)

// -------------------
// ReferenceSet
// -------------------

ReferenceSet::ReferenceSet()
    : DataSetBase("PacBio.DataSet.ReferenceSet", "ReferenceSet", XsdType::DATASETS)
{
}

// -------------------
// SubDataSets
// -------------------

SubDataSets::SubDataSets()
    : internal::DataSetListElement<DataSetBase>("DataSets", XsdType::DATASETS)
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

// -------------------
// TranscriptSet
// -------------------

TranscriptSet::TranscriptSet()
    : DataSetBase("PacBio.DataSet.TranscriptSet", "TranscriptSet", XsdType::DATASETS)
{
}

}  // namespace BAM
}  // namespace PacBio
