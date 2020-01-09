// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "RunMetadataParser.h"

#include <fstream>
#include <iostream>
#include <stdexcept>

#include <boost/algorithm/string.hpp>

#include <pbcopper/utility/StringUtils.h>

namespace PacBio {
namespace BAM {
namespace {

boost::optional<PacBio::BAM::AutomationParameters> AutomationParametersFromXml(
    const pugi::xml_node& paramsNode)
{
    if (!paramsNode) return boost::none;

    std::map<std::string, std::string> params;
    for (const auto& child : paramsNode.children()) {
        if (std::string{child.name()} == Element::AutomationParameter)
            params.emplace(child.attribute("Name").value(), child.attribute("SimpleValue").value());
    }
    return AutomationParameters{std::move(params)};
}

boost::optional<PacBio::BAM::BindingKit> BindingKitFromXml(const pugi::xml_node& kitNode)
{
    if (!kitNode) return boost::none;

    BindingKit kit;
    kit.PartNumber(kitNode.attribute(Element::PartNumber).value());
    return kit;
}

boost::optional<PacBio::BAM::ControlKit> ControlKitFromXml(const pugi::xml_node& kitNode)
{
    if (!kitNode) return boost::none;

    std::map<std::string, std::string> data;
    data.emplace(Element::PartNumber, kitNode.attribute(Element::PartNumber).value());

    const auto customSeqNode = kitNode.child(Element::CustomSequence);
    if (customSeqNode) {

        const auto lines = [](const std::string& input) {
            std::vector<std::string> result;
            size_t pos = 0;
            size_t found = input.find("\\n");
            while (found != std::string::npos) {
                result.push_back(input.substr(pos, found - pos));
                pos = found + 2;  // "\n"
                found = input.find("\\n", pos);
            }
            result.push_back(input.substr(pos));  // store last
            return result;
        }(std::string{customSeqNode.text().get()});

        if (lines.size() != 6) {
            throw std::runtime_error{
                "[pbbam] run metadata ERROR: malformatted CustomSequence node"};
        }
        data.emplace(Element::LeftAdapter, lines.at(1));
        data.emplace(Element::RightAdapter, lines.at(3));
        data.emplace(Element::Sequence, lines.at(5));
    }

    return ControlKit{std::move(data)};
}

boost::optional<PacBio::BAM::SequencingKitPlate> SequencingKitPlateFromXml(
    const pugi::xml_node& kitNode)
{
    if (!kitNode) return boost::none;

    SequencingKitPlate kit;
    kit.PartNumber(kitNode.attribute(Element::PartNumber).value());
    return kit;
}

boost::optional<PacBio::BAM::TemplatePrepKit> TemplatePrepKitFromXml(const pugi::xml_node& kitNode)
{
    if (!kitNode) return boost::none;

    std::map<std::string, std::string> data;
    data.emplace(Element::PartNumber, kitNode.attribute(Element::PartNumber).value());
    data.emplace(Element::LeftAdaptorSequence,
                 kitNode.child(Element::LeftAdaptorSequence).text().get());

    data.emplace(Element::LeftPrimerSequence,
                 kitNode.child(Element::LeftPrimerSequence).text().get());

    data.emplace(Element::RightAdaptorSequence,
                 kitNode.child(Element::RightAdaptorSequence).text().get());

    data.emplace(Element::RightPrimerSequence,
                 kitNode.child(Element::RightPrimerSequence).text().get());

    return TemplatePrepKit{std::move(data)};
}

CollectionMetadata SubreadSetCollection(const std::string& subreadSetName,
                                        const pugi::xml_node& subreadSetNode)
{
    const auto cmNode = subreadSetNode.child(Element::DataSetMetadata)
                            .child(Element::Collections)
                            .child(Element::CollectionMetadata);
    if (!cmNode)
        throw std::runtime_error{"[pbbam] run metadata ERROR: XML is missing expected elements"};

    return CollectionMetadata{
        subreadSetName, AutomationParametersFromXml(
                            cmNode.child(Element::Automation).child(Element::AutomationParameters)),
        BindingKitFromXml(cmNode.child(Element::BindingKit)),
        ControlKitFromXml(cmNode.child(Element::ControlKit)),
        SequencingKitPlateFromXml(cmNode.child(Element::SequencingKitPlate)),
        TemplatePrepKitFromXml(cmNode.child(Element::TemplatePrepKit))};
}

pugi::xml_node FetchSubreadSetsNode(const pugi::xml_document& doc)
{
    const auto rootNode = doc.document_element();
    if (!rootNode)
        throw std::runtime_error{"[pbbam] run metadata ERROR: could not fetch XML root node"};
    if (std::string{rootNode.name()} != Element::PacBioDataModel) {
        throw std::runtime_error{
            "[pbbam] run metadata ERROR: expected 'PacBioDataModel' as root node, instead found: " +
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
    if (loadResult.status != pugi::status_ok) {
        throw std::runtime_error{"[pbbam] run metadata ERROR: could not read XML, error code:" +
                                 std::to_string(loadResult.status)};
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