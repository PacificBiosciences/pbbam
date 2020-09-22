#include "PbbamInternalConfig.h"

#include "RunMetadataParser.h"

#include <fstream>
#include <istream>
#include <stdexcept>

#include <boost/algorithm/string.hpp>

#include <pbcopper/utility/StringUtils.h>

namespace PacBio {
namespace BAM {
namespace {

std::shared_ptr<internal::DataSetElement> MakeRunMetadataElement(const pugi::xml_node& xmlNode)
{
    std::string name = xmlNode.name();
    const auto foundColon = name.find(':');
    if (foundColon != std::string::npos) {
        name = name.substr(foundColon + 1);
    }

    const internal::FromInputXml fromInputXml;
    if (name == Element::Automation) return std::make_shared<Automation>(fromInputXml);
    if (name == Element::AutomationParameter)
        return std::make_shared<AutomationParameter>(fromInputXml);
    if (name == Element::AutomationParameters)
        return std::make_shared<AutomationParameters>(fromInputXml);
    if (name == Element::BindingKit) return std::make_shared<BindingKit>(fromInputXml);
    if (name == Element::Collections) return std::make_shared<Collections>(fromInputXml);
    if (name == Element::ControlKit) return std::make_shared<ControlKit>(fromInputXml);
    if (name == Element::SequencingKitPlate)
        return std::make_shared<SequencingKitPlate>(fromInputXml);
    if (name == Element::TemplatePrepKit) return std::make_shared<TemplatePrepKit>(fromInputXml);

    return std::make_shared<internal::DataSetElement>(name, internal::FromInputXml{});
}

void FromRunMetadataXml(const pugi::xml_node& xmlNode, internal::DataSetElement& parent)
{
    const std::string label = xmlNode.name();
    if (label.empty()) return;

    auto e = MakeRunMetadataElement(xmlNode);
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
        FromRunMetadataXml(childNode, *e.get());
    }

    parent.AddChild(e);
}

CollectionMetadata SubreadSetCollection(const std::string& subreadSetName,
                                        const pugi::xml_node& subreadSetNode)
{
    // find & initialize CollectionMetadata from node
    const auto cmNode = subreadSetNode.child(Element::DataSetMetadata)
                            .child(Element::Collections)
                            .child(Element::CollectionMetadata);
    if (!cmNode)
        throw std::runtime_error{"[pbbam] run metadata ERROR: XML is missing expected elements"};

    CollectionMetadata cm{subreadSetName};
    cm.Label(cmNode.name());

    // load element attributes
    auto attrIter = cmNode.attributes_begin();
    auto attrEnd = cmNode.attributes_end();
    for (; attrIter != attrEnd; ++attrIter) {
        const std::string name = attrIter->name();
        const std::string value = attrIter->value();
        cm.Attribute(name, value);
    }

    // load children, recursively
    auto childIter = cmNode.begin();
    auto childEnd = cmNode.end();
    for (; childIter != childEnd; ++childIter) {
        pugi::xml_node childNode = *childIter;
        FromRunMetadataXml(childNode, cm);
    }

    return cm;
}

pugi::xml_node FetchSubreadSetsNode(const pugi::xml_document& doc)
{
    const auto rootNode = doc.document_element();
    if (!rootNode)
        throw std::runtime_error{"[pbbam] run metadata ERROR: could not fetch XML root node"};
    if (std::string{rootNode.name()} != Element::PacBioDataModel) {
        throw std::runtime_error{
            "[pbbam] run metadata ERROR: expected 'PacBioDataModel' as root node, instead "
            "found: " +
            std::string{rootNode.name()}};
    }

    const auto result = rootNode.child(Element::ExperimentContainer)
                            .child(Element::Runs)
                            .child(Element::Run)
                            .child(Element::Outputs)
                            .child(Element::SubreadSets);
    if (!result)
        throw std::runtime_error{"[pbbam] run metadata ERROR: XML is missing expected elements"};
    return result;
}

std::map<std::string, CollectionMetadata> CollectionsFromXml(std::istream& in)
{
    std::map<std::string, CollectionMetadata> collections;

    pugi::xml_document doc;
    const pugi::xml_parse_result loadResult = doc.load(in);
    if (!loadResult) {
        throw std::runtime_error{
            "[pbbam] run metadata ERROR: could not read XML document\n"
            "  reason: " +
            std::string{loadResult.description()}};
    }

    const auto subreadSetsNode = FetchSubreadSetsNode(doc);
    for (const auto& subreadSetNode : subreadSetsNode) {
        const auto& subreadSetName = subreadSetNode.attribute("Name").value();
        collections.emplace(subreadSetName, SubreadSetCollection(subreadSetName, subreadSetNode));
    }

    return collections;
}

}  // namespace

CollectionMetadata RunMetadataParser::LoadCollection(const std::string& metadataXmlFn)
{
    std::ifstream in{metadataXmlFn};
    return LoadCollection(in);
}

CollectionMetadata RunMetadataParser::LoadCollection(std::istream& in)
{
    const auto& collections = CollectionsFromXml(in);
    // enforce only a single collection
    if (collections.size() != 1) {
        throw std::runtime_error{
            "[pbbam] run metadata ERROR: expected 1 SubreadSet, instead found: " +
            std::to_string(collections.size())};
    }
    return collections.begin()->second;
}

std::map<std::string, CollectionMetadata> RunMetadataParser::LoadCollections(
    const std::string& runMetadataXmlFn)
{
    std::ifstream in{runMetadataXmlFn};
    return LoadCollections(in);
}

std::map<std::string, CollectionMetadata> RunMetadataParser::LoadCollections(std::istream& in)
{
    return CollectionsFromXml(in);
}

}  // namespace BAM
}  // namespace PacBio
