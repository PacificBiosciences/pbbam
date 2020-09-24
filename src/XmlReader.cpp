#include "PbbamInternalConfig.h"

#include "XmlReader.h"

#include <cassert>

#include <istream>
#include <memory>
#include <stdexcept>
#include <typeinfo>
#include <vector>

#include <pbbam/StringUtilities.h>

#include "pugixml/pugixml.hpp"

using DataSetElement = PacBio::BAM::internal::DataSetElement;
using FromInputXml = PacBio::BAM::internal::FromInputXml;

namespace PacBio {
namespace BAM {
namespace {

std::unique_ptr<DataSetBase> MakeDataSetBase(const pugi::xml_node& xmlNode)
{
    const FromInputXml fromInputXml{};
    std::string name = xmlNode.name();
    const auto foundColon = name.find(':');
    if (foundColon != std::string::npos) {
        name = name.substr(foundColon + 1);
    }

    const auto type = ElementTypeFromName(name);
    switch (type) {
        case XmlElementType::ALIGNMENT_SET:
            return std::make_unique<AlignmentSet>(fromInputXml);
        case XmlElementType::BARCODE_SET:
            return std::make_unique<BarcodeSet>(fromInputXml);
        case XmlElementType::CONSENSUS_ALIGNMENT_SET:
            return std::make_unique<ConsensusAlignmentSet>(fromInputXml);
        case XmlElementType::CONSENSUS_READ_SET:
            return std::make_unique<ConsensusReadSet>(fromInputXml);
        case XmlElementType::CONTIG_SET:
            return std::make_unique<ContigSet>(fromInputXml);
        case XmlElementType::HDF_SUBREAD_SET:
            return std::make_unique<HdfSubreadSet>(fromInputXml);
        case XmlElementType::REFERENCE_SET:
            return std::make_unique<ReferenceSet>(fromInputXml);
        case XmlElementType::SUBREAD_SET:
            return std::make_unique<SubreadSet>(fromInputXml);
        case XmlElementType::TRANSCRIPT_SET:
            return std::make_unique<TranscriptSet>(fromInputXml);
        case XmlElementType::TRANSCRIPT_ALIGNMENT_SET:
            return std::make_unique<TranscriptAlignmentSet>(fromInputXml);
        case XmlElementType::GENERIC_DATASET:
            return std::make_unique<DataSetBase>(fromInputXml);
        default:
            // unreachable
            throw std::runtime_error{"[pbbam] XML reader ERROR: unknown data set label: " + name};
    }
}

std::shared_ptr<DataSetElement> MakeElement(const pugi::xml_node& xmlNode)
{
    std::string name = xmlNode.name();
    const auto foundColon = name.find(':');
    if (foundColon != std::string::npos) {
        name = name.substr(foundColon + 1);
    }

    const FromInputXml fromInputXml{};
    const auto type = ElementTypeFromName(name);
    switch (type) {
        case XmlElementType::AUTOMATION:
            return std::make_shared<Automation>(fromInputXml);
        case XmlElementType::AUTOMATION_PARAMETER:
            return std::make_shared<AutomationParameter>(fromInputXml);
        case XmlElementType::AUTOMATION_PARAMETERS:
            return std::make_shared<AutomationParameters>(fromInputXml);
        case XmlElementType::BINDING_KIT:
            return std::make_shared<BindingKit>(fromInputXml);
        case XmlElementType::BIOSAMPLE:
            return std::make_shared<BioSample>("", fromInputXml);
        case XmlElementType::BIOSAMPLES:
            return std::make_shared<BioSamples>(fromInputXml);
        case XmlElementType::COLLECTIONS:
            return std::make_shared<Collections>(fromInputXml);
        case XmlElementType::COLLECTION_METADATA:
            return std::make_shared<CollectionMetadata>(fromInputXml);
        case XmlElementType::CONTROL_KIT:
            return std::make_shared<ControlKit>(fromInputXml);
        case XmlElementType::DATASET_METADATA:
            return std::make_shared<DataSetMetadata>(fromInputXml);
        case XmlElementType::DNA_BARCODE:
            return std::make_shared<DNABarcode>("", fromInputXml);
        case XmlElementType::DNA_BARCODES:
            return std::make_shared<DNABarcodes>(fromInputXml);
        case XmlElementType::EXTENSION:
            return std::make_shared<ExtensionElement>(fromInputXml);
        case XmlElementType::EXTENSIONS:
            return std::make_shared<Extensions>(fromInputXml);
        case XmlElementType::EXTERNAL_RESOURCE:
            return std::make_shared<ExternalResource>("", "", fromInputXml);
        case XmlElementType::EXTERNAL_RESOURCES:
            return std::make_shared<ExternalResources>(fromInputXml);
        case XmlElementType::FILE_INDEX:
            return std::make_shared<FileIndex>("", "", fromInputXml);
        case XmlElementType::FILE_INDICES:
            return std::make_shared<FileIndices>(fromInputXml);
        case XmlElementType::FILTER:
            return std::make_shared<Filter>(fromInputXml);
        case XmlElementType::FILTERS:
            return std::make_shared<Filters>(fromInputXml);
        case XmlElementType::PARENT_TOOL:
            return std::make_shared<ParentTool>(fromInputXml);
        case XmlElementType::PPACONFIG:
            return std::make_shared<PPAConfig>(fromInputXml);
        case XmlElementType::PROPERTY:
            return std::make_shared<Property>("", "", "", fromInputXml);
        case XmlElementType::PROPERTIES:
            return std::make_shared<Properties>(fromInputXml);
        case XmlElementType::PROVENANCE:
            return std::make_shared<Provenance>(fromInputXml);
        case XmlElementType::SEQUENCING_KIT_PLATE:
            return std::make_shared<SequencingKitPlate>(fromInputXml);
        case XmlElementType::TEMPLATE_PREP_KIT:
            return std::make_shared<TemplatePrepKit>(fromInputXml);

        case XmlElementType::ALIGNMENT_SET:
            return std::make_shared<AlignmentSet>(fromInputXml);
        case XmlElementType::BARCODE_SET:
            return std::make_shared<BarcodeSet>(fromInputXml);
        case XmlElementType::CONSENSUS_ALIGNMENT_SET:
            return std::make_shared<ConsensusAlignmentSet>(fromInputXml);
        case XmlElementType::CONSENSUS_READ_SET:
            return std::make_shared<ConsensusReadSet>(fromInputXml);
        case XmlElementType::CONTIG_SET:
            return std::make_shared<ContigSet>(fromInputXml);
        case XmlElementType::HDF_SUBREAD_SET:
            return std::make_shared<HdfSubreadSet>(fromInputXml);
        case XmlElementType::SUBREAD_SET:
            return std::make_shared<SubreadSet>(fromInputXml);
        case XmlElementType::REFERENCE_SET:
            return std::make_shared<ReferenceSet>(fromInputXml);
        case XmlElementType::TRANSCRIPT_SET:
            return std::make_shared<TranscriptSet>(fromInputXml);
        case XmlElementType::TRANSCRIPT_ALIGNMENT_SET:
            return std::make_shared<TranscriptAlignmentSet>(fromInputXml);
        case XmlElementType::SUBDATASETS:
            return std::make_shared<SubDataSets>(fromInputXml);
        case XmlElementType::GENERIC_DATASET:
            return std::make_shared<DataSetBase>(fromInputXml);
        case XmlElementType::GENERIC_ELEMENT:
            return std::make_shared<DataSetElement>(name, fromInputXml);
        default:
            // unreachable
            throw std::runtime_error{"[pbbam] XML reader ERROR: unknown data element label: " +
                                     name};
    }
}

void UpdateRegistry(const std::string& attributeName, const std::string& attributeValue,
                    NamespaceRegistry& registry)
{
    std::vector<std::string> nameParts = Split(attributeName, ':');
    assert(!nameParts.empty());
    if (nameParts.size() > 2)
        throw std::runtime_error{"[pbbam] XML reader ERROR: malformed xmlns attribute: " +
                                 attributeName};

    const bool isDefault = (nameParts.size() == 1);
    const XsdType xsd = registry.XsdForUri(attributeValue);

    if (isDefault)
        registry.SetDefaultXsd(xsd);
    else {
        assert(nameParts.size() == 2);
        const std::string& name = nameParts.at(1);
        const std::string& uri = attributeValue;
        NamespaceInfo namespaceInfo(name, uri);
        registry.Register(xsd, namespaceInfo);
    }
}

void FromXml(const pugi::xml_node& xmlNode, DataSetElement& parent)
{
    // ignore non-named XML nodes
    //
    // pugi::xml separates XML parts into more node types than we use
    //
    const std::string label = xmlNode.name();
    if (label.empty()) return;

    auto e = MakeElement(xmlNode);
    e->Label(xmlNode.name());
    e->Text(xmlNode.text().get());

    // iterate attributes
    auto attrIter = xmlNode.attributes_begin();
    auto attrEnd = xmlNode.attributes_end();
    for (; attrIter != attrEnd; ++attrIter)
        e->Attribute(attrIter->name(), attrIter->value());

    // iterate children, recursively building up subtree
    auto childIter = xmlNode.begin();
    auto childEnd = xmlNode.end();
    for (; childIter != childEnd; ++childIter) {
        pugi::xml_node childNode = *childIter;
        FromXml(childNode, *e.get());
    }

    parent.AddChild(e);
}

}  // namespace

std::unique_ptr<DataSetBase> XmlReader::FromStream(std::istream& in)
{
    pugi::xml_document doc;
    const pugi::xml_parse_result loadResult = doc.load(in);
    if (!loadResult) {
        throw std::runtime_error{
            "[pbbam] XML reader ERROR: could not read XML document\n"
            "  reason: " +
            std::string{loadResult.description()}};
    }

    // parse top-level attributes
    pugi::xml_node rootNode = doc.document_element();
    if (rootNode == pugi::xml_node())
        throw std::runtime_error{"[pbbam] XML reader ERROR: could not fetch XML root node"};

    // create dataset matching type strings
    auto dataset = MakeDataSetBase(rootNode);
    dataset->Label(rootNode.name());

    // iterate attributes, capture namespace info
    const std::string xmlnsPrefix("xmlns");
    auto attrIter = rootNode.attributes_begin();
    auto attrEnd = rootNode.attributes_end();
    for (; attrIter != attrEnd; ++attrIter) {
        const std::string name = attrIter->name();
        const std::string value = attrIter->value();
        dataset->Attribute(name, value);
        if (name.find(xmlnsPrefix) == 0) {
            UpdateRegistry(name, value, dataset->Namespaces());
        }
    }

    // iterate children, recursively building up subtree
    auto childIter = rootNode.begin();
    auto childEnd = rootNode.end();
    for (; childIter != childEnd; ++childIter) {
        pugi::xml_node childNode = *childIter;
        FromXml(childNode, *dataset.get());
    }

    return dataset;
}

}  // namespace BAM
}  // namespace PacBio
