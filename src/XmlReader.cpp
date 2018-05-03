// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "XmlReader.h"

#include <cassert>
#include <fstream>
#include <iostream>
#include <memory>
#include <vector>

#include "StringUtils.h"
#include "pugixml/pugixml.hpp"

namespace PacBio {
namespace BAM {
namespace internal {

static void UpdateRegistry(const std::string& attributeName, const std::string& attributeValue,
                           NamespaceRegistry& registry)
{
    std::vector<std::string> nameParts = Split(attributeName, ':');
    assert(!nameParts.empty());
    if (nameParts.size() > 2)
        throw std::runtime_error{"malformed xmlns attribute: " + attributeName};

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

static void FromXml(const pugi::xml_node& xmlNode, DataSetElement& parent)
{
    // ignore non-named XML nodes
    //
    // pugi::xml separates XML parts into more node types than we use
    //
    const std::string label = xmlNode.name();
    if (label.empty()) return;

    // label & text
    DataSetElement e(xmlNode.name(), FromInputXml());
    e.Text(xmlNode.text().get());

    // iterate attributes
    auto attrIter = xmlNode.attributes_begin();
    auto attrEnd = xmlNode.attributes_end();
    for (; attrIter != attrEnd; ++attrIter)
        e.Attribute(attrIter->name(), attrIter->value());

    // iterate children, recursively building up subtree
    auto childIter = xmlNode.begin();
    auto childEnd = xmlNode.end();
    for (; childIter != childEnd; ++childIter) {
        pugi::xml_node childNode = *childIter;
        FromXml(childNode, e);
    }

    // add our element to its parent
    parent.AddChild(e);
}

std::unique_ptr<DataSetBase> XmlReader::FromStream(std::istream& in)
{
    pugi::xml_document doc;
    const pugi::xml_parse_result loadResult = doc.load(in);
    if (loadResult.status != pugi::status_ok)
        throw std::runtime_error{"could not read XML file, error code:" +
                                 std::to_string(loadResult.status)};

    // parse top-level attributes
    pugi::xml_node rootNode = doc.document_element();
    if (rootNode == pugi::xml_node()) throw std::runtime_error{"could not fetch XML root node"};

    // create dataset matching type strings
    std::unique_ptr<DataSetBase> dataset(new DataSetBase);
    dataset->Label(rootNode.name());

    // iterate attributes, capture namespace info
    const std::string xmlnsPrefix("xmlns");
    auto attrIter = rootNode.attributes_begin();
    auto attrEnd = rootNode.attributes_end();
    for (; attrIter != attrEnd; ++attrIter) {
        const std::string name = attrIter->name();
        const std::string value = attrIter->value();
        dataset->Attribute(name, value);

        if (name.find(xmlnsPrefix) == 0) UpdateRegistry(name, value, dataset->Namespaces());
    }

    // iterate children, recursively building up subtree
    auto childIter = rootNode.begin();
    auto childEnd = rootNode.end();
    for (; childIter != childEnd; ++childIter) {
        pugi::xml_node childNode = *childIter;
        internal::FromXml(childNode, *dataset.get());
    }

    return dataset;
}

}  // namespace internal
}  // namespace BAM
}  // namespace PacBio
