// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "XmlWriter.h"

#include <cstddef>

#include <map>
#include <ostream>
#include <stdexcept>

#include "pbbam/DataSet.h"

#include "FileUtils.h"
#include "pugixml/pugixml.hpp"

using DataSetElement = PacBio::BAM::internal::DataSetElement;

namespace PacBio {
namespace BAM {
namespace {

std::string Prefix(const std::string& input)
{
    const auto colonFound = input.find(':');
    if (colonFound == std::string::npos || colonFound == 0) return std::string();
    return input.substr(0, colonFound);
}

std::string OutputName(const DataSetElement& node, const NamespaceRegistry& registry)
{
    // if from input XML, respect the namespaces given
    if (node.IsVerbatimLabel())
        return node.QualifiedNameLabel();

    else if (node.LocalNameLabel().to_string() == "Collections")
        return "Collections";

    // otherwise, probably user-generated
    else {
        // if no namespace prefix, prepend the appropriate one & return
        if (node.PrefixLabel().empty()) {
            static const std::string colon = ":";
            auto xsdType = node.Xsd();
            if (xsdType == XsdType::NONE)
                xsdType = registry.XsdForElement(node.LocalNameLabel().to_string());
            return registry.Namespace(xsdType).Name() + colon + node.LocalNameLabel().to_string();
        }
        // otherwise, has prefix - return full name
        else
            return node.QualifiedNameLabel();
    }
}

// (admitted hack to) ensure order of CollectionMetadata children
void CollectionMetadataToXml(const DataSetElement& node, const NamespaceRegistry& registry,
                             std::map<XsdType, std::string>& xsdPrefixesUsed,
                             pugi::xml_node& parentXml, const DataSetBase& dataset,
                             DataSetPathMode pathMode);

void ToXml(const DataSetElement& node, const NamespaceRegistry& registry,
           std::map<XsdType, std::string>& xsdPrefixesUsed, pugi::xml_node& parentXml,
           const DataSetBase& dataset, DataSetPathMode pathMode)
{
    // create child of parent, w/ label & text
    const auto label = OutputName(node, registry);
    if (label.empty()) return;  // error?
    auto xmlNode = parentXml.append_child(label.c_str());

    if (!node.Text().empty()) xmlNode.text().set(node.Text().c_str());

    // store XSD type for later
    const auto prefix = Prefix(label);
    if (!prefix.empty()) {
        if (label.find("pbmeta:AutomationParameter") == std::string::npos) {
            xsdPrefixesUsed[node.Xsd()] = prefix;
        }
    }

    // add attributes
    const bool resolveFilePaths = (pathMode == DataSetPathMode::ABSOLUTE);
    for (const auto& attribute : node.Attributes()) {
        const auto& name = attribute.first;
        if (name.empty()) continue;

        auto attr = xmlNode.append_attribute(name.c_str());
        std::string value = attribute.second.c_str();
        // "absolutize" any paths, except relative paths from verbatim input XML
        if (!dataset.FromInputXml() && resolveFilePaths && name == "ResourceId")
            value = FileUtils::ResolvedFilePath(value, dataset.Path());
        attr.set_value(value.c_str());
    }

    // additional stuff later? (e.g. comments)

    // iterate children, recursively building up subtree
    for (const auto& child : node.Children()) {
        // ensure order CollectionMetadata children
        if (child->QualifiedNameLabel().find("CollectionMetadata") != std::string::npos)
            CollectionMetadataToXml(*child, registry, xsdPrefixesUsed, xmlNode, dataset, pathMode);
        else
            ToXml(*child, registry, xsdPrefixesUsed, xmlNode, dataset, pathMode);
    }
}

// (admitted hack to) ensure order of CollectionMetadata children
void CollectionMetadataToXml(const DataSetElement& node, const NamespaceRegistry& registry,
                             std::map<XsdType, std::string>& xsdPrefixesUsed,
                             pugi::xml_node& parentXml, const DataSetBase& dataset,
                             DataSetPathMode pathMode)
{
    // create child of parent, w/ label & text
    const auto label = OutputName(node, registry);
    if (label.empty()) return;  // error?
    auto xmlNode = parentXml.append_child(label.c_str());

    if (!node.Text().empty()) xmlNode.text().set(node.Text().c_str());

    // store XSD type for later
    const auto prefix = Prefix(label);
    if (!prefix.empty()) {
        if (label.find("pbmeta:AutomationParameter") == std::string::npos) {
            xsdPrefixesUsed[node.Xsd()] = prefix;
        }
    }

    // add attributes
    const bool resolveFilePaths = (pathMode == DataSetPathMode::ABSOLUTE);
    for (const auto& attribute : node.Attributes()) {
        const auto& name = attribute.first;
        if (name.empty()) continue;

        auto attr = xmlNode.append_attribute(name.c_str());
        std::string value = attribute.second.c_str();
        // "absolutize" any paths, except relative paths from verbatim input XML
        if (!dataset.FromInputXml() && resolveFilePaths && name == "ResourceId")
            value = FileUtils::ResolvedFilePath(value, dataset.Path());
        attr.set_value(value.c_str());
    }

    // force child ordering
    static const std::vector<std::string> childOrder{
        "MultiJobId",         "ConsensusReadSetRef", "InstCtrlVer",
        "SigProcVer",         "RunDetails",          "Movie",
        "WellSample",         "Automation",          "CollectionNumber",
        "CellIndex",          "SetNumber",           "CellPac",
        "ControlKit",         "TemplatePrepKit",     "BindingKit",
        "SequencingKitPlate", "WashKitPlate",        "Primary",
        "PPAConfig",          "Secondary",           "UserDefinedFields",
        "ExpirationData",     "ComponentVersions"};
    std::set<std::string> seen;
    for (const auto& name : childOrder) {
        if (node.HasChild(name)) {
            seen.insert(node.QualifiedNameLabel());
            const DataSetElement& child = node.Child<DataSetElement>(name);
            ToXml(child, registry, xsdPrefixesUsed, xmlNode, dataset, pathMode);
        }
    }

    // any remaining nodes will be unordered w.r.t. children above,
    // but at least they're not lost
    for (const auto& child : node.Children()) {
        const auto insertResult = seen.insert(child->QualifiedNameLabel());
        const auto alreadySeen = insertResult.second;
        if (!alreadySeen) {
            ToXml(*child, registry, xsdPrefixesUsed, xmlNode, dataset, pathMode);
        }
    }

    return;
}

}  // namespace

void XmlWriter::ToStream(const DataSetBase& dataset, std::ostream& out, DataSetPathMode pathMode)
{
    pugi::xml_document doc;

    const auto& registry = dataset.Namespaces();

    // create top-level dataset XML node
    const auto label = OutputName(dataset, registry);
    if (label.empty())
        throw std::runtime_error{"[pbbam] XML writer ERROR: could not convert dataset node to XML"};
    auto root = doc.append_child(label.c_str());

    const auto& text = dataset.Text();
    if (!text.empty()) root.text().set(text.c_str());

    // add top-level attributes
    for (const auto& attribute : dataset.Attributes()) {
        const auto& name = attribute.first;
        const auto& value = attribute.second;
        if (name.empty()) continue;
        auto attr = root.append_attribute(name.c_str());
        attr.set_value(value.c_str());
    }

    std::map<XsdType, std::string> xsdPrefixesUsed;
    xsdPrefixesUsed[dataset.Xsd()] = Prefix(label);
    if (dataset.Metadata().HasChild("Collections")) {
        xsdPrefixesUsed[XsdType::COLLECTION_METADATA] = "pbmeta";
    }

    // iterate children, recursively building up subtree
    for (const auto& child : dataset.Children())
        ToXml(*child, registry, xsdPrefixesUsed, root, dataset, pathMode);

    // write XML to stream
    auto decl = doc.prepend_child(pugi::node_declaration);
    decl.append_attribute("version") = "1.0";
    decl.append_attribute("encoding") = "utf-8";

    // add XSD namespace attributes
    auto xmlnsDefaultAttribute = root.attribute("xmlns");
    if (xmlnsDefaultAttribute.empty()) {
        xmlnsDefaultAttribute = root.append_attribute("xmlns");
        xmlnsDefaultAttribute.set_value(registry.DefaultNamespace().Uri().c_str());
    }
    auto xsiAttribute = root.attribute("xmlns:xsi");
    if (xsiAttribute.empty()) {
        xsiAttribute = root.append_attribute("xmlns:xsi");
        xsiAttribute.set_value("http://www.w3.org/2001/XMLSchema-instance");
    }
    auto xsiSchemaLocationAttribute = root.attribute("xsi:schemaLocation");
    if (xsiSchemaLocationAttribute.empty()) {
        xsiSchemaLocationAttribute = root.append_attribute("xsi:schemaLocation");
        xsiSchemaLocationAttribute.set_value(registry.DefaultNamespace().Uri().c_str());
    }

    static const std::string xmlnsPrefix = "xmlns:";
    for (const auto& prefixIter : xsdPrefixesUsed) {
        const auto& xsdType = prefixIter.first;
        const auto& prefix = prefixIter.second;
        if (xsdType == XsdType::NONE || prefix.empty()) continue;

        const auto& nsInfo = registry.Namespace(xsdType);
        assert(nsInfo.Name() == prefix);
        const auto xmlnsName = xmlnsPrefix + prefix;
        auto xmlnsAttribute = root.attribute(xmlnsName.c_str());
        if (xmlnsAttribute.empty()) {
            xmlnsAttribute = root.append_attribute(xmlnsName.c_str());
            xmlnsAttribute.set_value(nsInfo.Uri().c_str());
        }
    }

    // "no escapes" to allow explicit ">" "<" comparison operators in filter parameters
    // we may remove this if/when comparison is separated from the value
    doc.save(out, "\t", pugi::format_default, pugi::encoding_utf8);
}

void XmlWriter::ToStream(const std::unique_ptr<DataSetBase>& dataset, std::ostream& out,
                         DataSetPathMode pathMode)
{
    ToStream(*dataset.get(), out, pathMode);
}

}  // namespace BAM
}  // namespace PacBio
