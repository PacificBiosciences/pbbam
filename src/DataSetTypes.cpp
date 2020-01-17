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
    {"BioSample",              ElementType::BIOSAMPLE},
    {"BioSamples",             ElementType::BIOSAMPLES},
    {"Collections",            ElementType::COLLECTIONS},
    {"CollectionMetadata",     ElementType::COLLECTION_METADATA},
    {"DNABarcode",             ElementType::DNA_BARCODE},
    {"DNABarcodes",            ElementType::DNA_BARCODES},
    {"ExtensionElement",       ElementType::EXTENSION},
    {"Extensions",             ElementType::EXTENSIONS},
    {"ExternalResource",       ElementType::EXTERNAL_RESOURCE},
    {"ExternalResources",      ElementType::EXTERNAL_RESOURCES},
    {"FileIndex",              ElementType::FILE_INDEX},
    {"FileIndices",            ElementType::FILE_INDICES},
    {"Filter",                 ElementType::FILTER},
    {"Filters",                ElementType::FILTERS},
    {"ParentTool",             ElementType::PARENT_TOOL},
    {"PPAConfig",              ElementType::PPACONFIG},
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

}  // namespace

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

// -------------------
// BioSample
// -------------------

BioSample::BioSample(const std::string& name) : DataSetElement("BioSample", XsdType::SAMPLE_INFO)
{
    Name(name);
}

BioSample::BioSample(const std::string& name, const internal::FromInputXml& fromInputXml)
    : DataSetElement("", fromInputXml, XsdType::SAMPLE_INFO)
{
    Name(name);
}

DEFINE_ACCESSORS(BioSample, DNABarcodes, DNABarcodes)

BioSample& BioSample::DNABarcodes(const PacBio::BAM::DNABarcodes& barcodes)
{
    DNABarcodes() = barcodes;
    return *this;
}

const std::string& BioSample::Name() const { return Attribute("Name"); }

std::string& BioSample::Name() { return Attribute("Name"); }

BioSample& BioSample::Name(const std::string& name)
{
    Attribute("Name", name);
    return *this;
}

// -------------------
// BioSamples
// -------------------

BioSamples::BioSamples() : DataSetElement("BioSamples", XsdType::SAMPLE_INFO) {}

BioSamples::BioSamples(const internal::FromInputXml& fromInputXml)
    : DataSetElement("", fromInputXml, XsdType::SAMPLE_INFO)
{
}

void BioSamples::Add(const BioSample& sample) { AddChild(sample); }

void BioSamples::Remove(const BioSample& sample) { RemoveChild(sample); }

BioSamples::iterator_type BioSamples::begin() { return BioSamples::iterator_type(this, 0); }

BioSamples::const_iterator_type BioSamples::begin() const { return cbegin(); }

BioSamples::const_iterator_type BioSamples::cbegin() const
{
    return BioSamples::const_iterator_type(this, 0);
}

BioSamples::iterator_type BioSamples::end()
{
    return BioSamples::iterator_type(this, NumChildren());
}

BioSamples::const_iterator_type BioSamples::end() const { return cend(); }

BioSamples::const_iterator_type BioSamples::cend() const
{
    return BioSamples::const_iterator_type(this, NumChildren());
}

const BioSamples::value_type& BioSamples::operator[](size_t index) const
{
    return dynamic_cast<const BioSamples::value_type&>(*(children_.at(index).get()));
}

BioSamples::value_type& BioSamples::operator[](size_t index)
{
    return dynamic_cast<BioSamples::value_type&>(*(children_.at(index).get()));
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
    , path_(FileUtils::CurrentWorkingDirectory())
{
}

DataSetBase::DataSetBase(const internal::FromInputXml& fromInputXml)
    : StrictEntityType("", "DataSet", fromInputXml, XsdType::DATASETS)
    , path_(FileUtils::CurrentWorkingDirectory())
{
}

DataSetBase::DataSetBase(const std::string& metatype, const std::string& label, const XsdType& xsd)
    : StrictEntityType(metatype, label, xsd), path_(FileUtils::CurrentWorkingDirectory())
{
}

DataSetBase::DataSetBase(const std::string& metatype, const std::string& label,
                         const internal::FromInputXml& fromInputXml, const XsdType& xsd)
    : StrictEntityType(metatype, label, fromInputXml, xsd)
    , path_(FileUtils::CurrentWorkingDirectory())
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

bool DataSetBase::FromInputXml() const { return fromInputXml_; }

void DataSetBase::FromInputXml(bool ok) { fromInputXml_ = ok; }

DEFINE_ACCESSORS(DataSetBase, DataSetMetadata, Metadata)

DataSetBase& DataSetBase::Metadata(const PacBio::BAM::DataSetMetadata& metadata)
{
    Metadata() = metadata;
    return *this;
}

const NamespaceRegistry& DataSetBase::Namespaces() const { return registry_; }

NamespaceRegistry& DataSetBase::Namespaces() { return registry_; }

void DataSetBase::Path(const std::string& path) { path_ = path; }

const std::string& DataSetBase::Path() const { return path_; }

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
    result->path_ = path_;
    return result;
}

DataSetBase& DataSetBase::operator+=(const DataSetBase& other)
{
    // must be same dataset types (or 'other' must be generic)
    if (other.LocalNameLabel() != LocalNameLabel() && other.LocalNameLabel() != "DataSet")
        throw std::runtime_error{"[pbbam] dataset ERROR: cannot merge different dataset types"};

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
    throw std::runtime_error{"[pbbam] dataset ERROR: unsupported type: " + typeName};
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
    throw std::runtime_error{"[pbbam] dataset ERROR: unsupported type: " + typeName};
}

void DataSetBase::Save(const std::string& outputFilename, DataSetPathMode pathMode)
{
    DataSetIO::ToFile(*this, outputFilename, pathMode);
}

void DataSetBase::SaveToStream(std::ostream& out, DataSetPathMode pathMode)
{
    DataSetIO::ToStream(*this, out, pathMode);
}

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

DEFINE_ACCESSORS(DataSetMetadata, BioSamples, BioSamples)

DataSetMetadata& DataSetMetadata::BioSamples(const PacBio::BAM::BioSamples& samples)
{
    BioSamples() = samples;
    return *this;
}

const PacBio::BAM::CollectionMetadata& DataSetMetadata::CollectionMetadata() const
{
    const PacBio::BAM::Collections& collections = Child<PacBio::BAM::Collections>("Collections");
    assert(collections.Size() >= 1);
    const PacBio::BAM::CollectionMetadata& cm =
        collections.Child<PacBio::BAM::CollectionMetadata>(0);
    return cm;
}

PacBio::BAM::CollectionMetadata& DataSetMetadata::CollectionMetadata()
{
    PacBio::BAM::Collections& collections = Child<PacBio::BAM::Collections>("Collections");
    if (collections.Size() == 0) collections.AddChild(PacBio::BAM::CollectionMetadata{});
    PacBio::BAM::CollectionMetadata& cm = collections.Child<PacBio::BAM::CollectionMetadata>(0);
    return cm;
}

DataSetMetadata& DataSetMetadata::CollectionMetadata(
    const PacBio::BAM::CollectionMetadata& metadata)
{
    CollectionMetadata() = metadata;
    return *this;
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

const std::string& DataSetMetadata::NumRecords() const { return ChildText("NumRecords"); }

std::string& DataSetMetadata::NumRecords() { return ChildText("NumRecords"); }

DataSetMetadata& DataSetMetadata::NumRecords(const std::string& numRecords)
{
    ChildText("NumRecords", numRecords);
    return *this;
}

const std::string& DataSetMetadata::TotalLength() const { return ChildText("TotalLength"); }

std::string& DataSetMetadata::TotalLength() { return ChildText("TotalLength"); }

DataSetMetadata& DataSetMetadata::TotalLength(const std::string& totalLength)
{
    ChildText("TotalLength", totalLength);
    return *this;
}

// -------------------
// DNABarcode
// -------------------

DNABarcode::DNABarcode(const std::string& name) : DataSetElement("DNABarcode", XsdType::SAMPLE_INFO)
{
    Name(name);
    UniqueId(internal::GenerateUuid());
}

DNABarcode::DNABarcode(const std::string& name, const std::string& uuid)
    : DataSetElement("DNABarcode", XsdType::SAMPLE_INFO)
{
    Name(name);
    UniqueId(uuid);
}

DNABarcode::DNABarcode(const std::string& name, const internal::FromInputXml& fromInputXml)
    : DataSetElement("", fromInputXml, XsdType::SAMPLE_INFO)
{
    Name(name);
    UniqueId(internal::GenerateUuid());
}

DNABarcode::DNABarcode(const std::string& name, const std::string& uuid,
                       const internal::FromInputXml& fromInputXml)
    : DataSetElement("", fromInputXml, XsdType::SAMPLE_INFO)
{
    Name(name);
    UniqueId(uuid);
}

const std::string& DNABarcode::Name() const { return Attribute("Name"); }

std::string& DNABarcode::Name() { return Attribute("Name"); }

DNABarcode& DNABarcode::Name(const std::string& name)
{
    Attribute("Name", name);
    return *this;
}

const std::string& DNABarcode::UniqueId() const { return Attribute("UniqueId"); }

std::string& DNABarcode::UniqueId() { return Attribute("UniqueId"); }

DNABarcode& DNABarcode::UniqueId(const std::string& uuid)
{
    Attribute("UniqueId", uuid);
    return *this;
}

// -------------------
// DNABarcodes
// -------------------

DNABarcodes::DNABarcodes() : DataSetElement("DNABarcodes", XsdType::SAMPLE_INFO) {}

DNABarcodes::DNABarcodes(const internal::FromInputXml& fromInputXml)
    : DataSetElement("", fromInputXml, XsdType::SAMPLE_INFO)
{
}

void DNABarcodes::Add(const DNABarcode& barcode) { AddChild(barcode); }

void DNABarcodes::Remove(const DNABarcode& barcode) { RemoveChild(barcode); }

DNABarcodes::iterator_type DNABarcodes::begin() { return DNABarcodes::iterator_type(this, 0); }

DNABarcodes::const_iterator_type DNABarcodes::begin() const { return cbegin(); }

DNABarcodes::const_iterator_type DNABarcodes::cbegin() const
{
    return DNABarcodes::const_iterator_type(this, 0);
}

DNABarcodes::iterator_type DNABarcodes::end()
{
    return DNABarcodes::iterator_type(this, NumChildren());
}

DNABarcodes::const_iterator_type DNABarcodes::end() const { return cend(); }

DNABarcodes::const_iterator_type DNABarcodes::cend() const
{
    return DNABarcodes::const_iterator_type(this, NumChildren());
}

const DNABarcodes::value_type& DNABarcodes::operator[](size_t index) const
{
    return dynamic_cast<const DNABarcodes::value_type&>(*(children_.at(index).get()));
}

DNABarcodes::value_type& DNABarcodes::operator[](size_t index)
{
    return dynamic_cast<DNABarcodes::value_type&>(*(children_.at(index).get()));
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

Extensions::iterator_type Extensions::begin() { return Extensions::iterator_type(this, 0); }

Extensions::const_iterator_type Extensions::begin() const { return cbegin(); }

Extensions::const_iterator_type Extensions::cbegin() const
{
    return Extensions::const_iterator_type(this, 0);
}

Extensions::iterator_type Extensions::end()
{
    return Extensions::iterator_type(this, NumChildren());
}

Extensions::const_iterator_type Extensions::end() const { return cend(); }

Extensions::const_iterator_type Extensions::cend() const
{
    return Extensions::const_iterator_type(this, NumChildren());
}

const Extensions::value_type& Extensions::operator[](size_t index) const
{
    return dynamic_cast<const Extensions::value_type&>(*(children_.at(index).get()));
}

Extensions::value_type& Extensions::operator[](size_t index)
{
    return dynamic_cast<Extensions::value_type&>(*(children_.at(index).get()));
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

ExternalResources::iterator_type ExternalResources::begin()
{
    return ExternalResources::iterator_type(this, 0);
}

ExternalResources::const_iterator_type ExternalResources::begin() const { return cbegin(); }

ExternalResources::const_iterator_type ExternalResources::cbegin() const
{
    return ExternalResources::const_iterator_type(this, 0);
}

ExternalResources::iterator_type ExternalResources::end()
{
    return ExternalResources::iterator_type(this, NumChildren());
}

ExternalResources::const_iterator_type ExternalResources::end() const { return cend(); }

ExternalResources::const_iterator_type ExternalResources::cend() const
{
    return ExternalResources::const_iterator_type(this, NumChildren());
}

const ExternalResources::value_type& ExternalResources::operator[](size_t index) const
{
    return dynamic_cast<const ExternalResources::value_type&>(*(children_.at(index).get()));
}

ExternalResources::value_type& ExternalResources::operator[](size_t index)
{
    return dynamic_cast<ExternalResources::value_type&>(*(children_.at(index).get()));
}

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

FileIndices::iterator_type FileIndices::begin() { return FileIndices::iterator_type(this, 0); }

FileIndices::const_iterator_type FileIndices::begin() const { return cbegin(); }

FileIndices::const_iterator_type FileIndices::cbegin() const
{
    return FileIndices::const_iterator_type(this, 0);
}

FileIndices::iterator_type FileIndices::end()
{
    return FileIndices::iterator_type(this, NumChildren());
}

FileIndices::const_iterator_type FileIndices::end() const { return cend(); }

FileIndices::const_iterator_type FileIndices::cend() const
{
    return FileIndices::const_iterator_type(this, NumChildren());
}

const FileIndices::value_type& FileIndices::operator[](size_t index) const
{
    return dynamic_cast<const FileIndices::value_type&>(*(children_.at(index).get()));
}

FileIndices::value_type& FileIndices::operator[](size_t index)
{
    return dynamic_cast<FileIndices::value_type&>(*(children_.at(index).get()));
}

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

Filters::iterator_type Filters::begin() { return Filters::iterator_type(this, 0); }

Filters::const_iterator_type Filters::begin() const { return cbegin(); }

Filters::const_iterator_type Filters::cbegin() const
{
    return Filters::const_iterator_type(this, 0);
}

Filters::iterator_type Filters::end() { return Filters::iterator_type(this, NumChildren()); }

Filters::const_iterator_type Filters::end() const { return cend(); }

Filters::const_iterator_type Filters::cend() const
{
    return Filters::const_iterator_type(this, NumChildren());
}

const Filters::value_type& Filters::operator[](size_t index) const
{
    return dynamic_cast<const Filters::value_type&>(*(children_.at(index).get()));
}

Filters::value_type& Filters::operator[](size_t index)
{
    return dynamic_cast<Filters::value_type&>(*(children_.at(index).get()));
}

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

Properties::iterator_type Properties::begin() { return Properties::iterator_type(this, 0); }

Properties::const_iterator_type Properties::begin() const { return cbegin(); }

Properties::const_iterator_type Properties::cbegin() const
{
    return Properties::const_iterator_type(this, 0);
}

Properties::iterator_type Properties::end()
{
    return Properties::iterator_type(this, NumChildren());
}

Properties::const_iterator_type Properties::end() const { return cend(); }

Properties::const_iterator_type Properties::cend() const
{
    return Properties::const_iterator_type(this, NumChildren());
}

const Properties::value_type& Properties::operator[](size_t index) const
{
    return dynamic_cast<const Properties::value_type&>(*(children_.at(index).get()));
}

Properties::value_type& Properties::operator[](size_t index)
{
    return dynamic_cast<Properties::value_type&>(*(children_.at(index).get()));
}

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

const std::string& Property::Name() const { return Attribute("Name"); }

std::string& Property::Name() { return Attribute("Name"); }

Property& Property::Name(const std::string& name)
{
    Attribute("Name", name);
    return *this;
}

const std::string& Property::Operator() const { return Attribute("Operator"); }

std::string& Property::Operator() { return Attribute("Operator"); }

Property& Property::Operator(const std::string& op)
{
    Attribute("Operator", op);
    return *this;
}

const std::string& Property::Value() const { return Attribute("Value"); }

std::string& Property::Value() { return Attribute("Value"); }

Property& Property::Value(const std::string& value)
{
    Attribute("Value", value);
    return *this;
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

const std::string& Provenance::CreatedBy() const { return Attribute("CreatedBy"); }

std::string& Provenance::CreatedBy() { return Attribute("CreatedBy"); }

Provenance& Provenance::CreatedBy(const std::string& createdBy)
{
    Attribute("CreatedBy", createdBy);
    return *this;
}

const std::string& Provenance::CommonServicesInstanceId() const
{
    return ChildText("CommonServicesInstanceId");
}

std::string& Provenance::CommonServicesInstanceId()
{
    return ChildText("CommonServicesInstanceId");
}

Provenance& Provenance::CommonServicesInstanceId(const std::string& id)
{
    ChildText("CommonServicesInstanceId", id);
    return *this;
}

const std::string& Provenance::CreatorUserId() const { return ChildText("CreatorUserId"); }

std::string& Provenance::CreatorUserId() { return ChildText("CreatorUserId"); }

Provenance& Provenance::CreatorUserId(const std::string& id)
{
    ChildText("CreatorUserId", id);
    return *this;
}

const std::string& Provenance::ParentJobId() const { return ChildText("ParentJobId"); }

std::string& Provenance::ParentJobId() { return ChildText("ParentJobId"); }

Provenance& Provenance::ParentJobId(const std::string& id)
{
    ChildText("ParentJobId", id);
    return *this;
}

Provenance& Provenance::ParentTool(const PacBio::BAM::ParentTool& tool)
{
    ParentTool() = tool;
    return *this;
}

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

SubDataSets::iterator_type SubDataSets::begin() { return SubDataSets::iterator_type(this, 0); }

SubDataSets::const_iterator_type SubDataSets::begin() const { return cbegin(); }

SubDataSets::const_iterator_type SubDataSets::cbegin() const
{
    return SubDataSets::const_iterator_type(this, 0);
}

SubDataSets::iterator_type SubDataSets::end()
{
    return SubDataSets::iterator_type(this, NumChildren());
}

SubDataSets::const_iterator_type SubDataSets::end() const { return cend(); }

SubDataSets::const_iterator_type SubDataSets::cend() const
{
    return SubDataSets::const_iterator_type(this, NumChildren());
}

const SubDataSets::value_type& SubDataSets::operator[](size_t index) const
{
    return dynamic_cast<const SubDataSets::value_type&>(*(children_.at(index).get()));
}

SubDataSets::value_type& SubDataSets::operator[](size_t index)
{
    return dynamic_cast<SubDataSets::value_type&>(*(children_.at(index).get()));
}

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
