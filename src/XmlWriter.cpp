// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "XmlWriter.h"
#include "pbbam/DataSet.h"
#include "pugixml/pugixml.hpp"
#include <cstddef>
#include <fstream>
#include <iostream>
#include <map>

namespace PacBio {
namespace BAM {
namespace internal {

static
std::string Prefix(const std::string& input)
{
    const auto colonFound = input.find(':');
    if (colonFound == std::string::npos || colonFound == 0)
        return std::string();
    return input.substr(0, colonFound);
}

static
std::string OutputName(const DataSetElement& node,
                       const NamespaceRegistry& registry)
{
    // if from input XML, respect the namespaces given
    if (node.IsVerbatimLabel()) 
        return node.QualifiedNameLabel();

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

static
void ToXml(const DataSetElement& node,
           const NamespaceRegistry& registry,
           std::map<XsdType, std::string>& xsdPrefixesUsed,
           pugi::xml_node& parentXml)
{
    // create child of parent, w/ label & text
    const auto label = OutputName(node, registry);
    if (label.empty())
        return; // error?
    auto xmlNode = parentXml.append_child(label.c_str());

    if (!node.Text().empty())
        xmlNode.text().set(node.Text().c_str());

    // store XSD type for later
    const auto prefix = Prefix(label);
    if (!prefix.empty())
        xsdPrefixesUsed[node.Xsd()] = prefix;

    // add attributes
    for (const auto& attribute : node.Attributes())
    {
        const auto& name = attribute.first;
        if (name.empty())
            continue;
        auto attr = xmlNode.append_attribute(name.c_str());
        attr.set_value(attribute.second.c_str());
    }

    // additional stuff later? (e.g. comments)

    // iterate children, recursively building up subtree
    for (const auto& child : node.Children())
        ToXml(child, registry, xsdPrefixesUsed, xmlNode);
}

void XmlWriter::ToStream(const DataSetBase& dataset,
                         std::ostream& out)
{
    pugi::xml_document doc;

    const auto& registry = dataset.Namespaces();

    // create top-level dataset XML node
    const auto label = internal::OutputName(dataset, registry);
    if (label.empty())
        throw std::runtime_error("could not convert dataset node to XML");
    auto root = doc.append_child(label.c_str());

    const auto& text = dataset.Text();
    if (!text.empty())
        root.text().set(text.c_str());

    // add top-level attributes
    for (const auto& attribute : dataset.Attributes()) {
        const auto& name = attribute.first;
        const auto& value = attribute.second;
        if (name.empty())
            continue;
        auto attr = root.append_attribute(name.c_str());
        attr.set_value(value.c_str());
    }

    std::map<XsdType, std::string> xsdPrefixesUsed;
    xsdPrefixesUsed[dataset.Xsd()] = Prefix(label);

    // iterate children, recursively building up subtree
    for (const auto& child : dataset.Children())
        ToXml(child, registry, xsdPrefixesUsed, root);

    // write XML to stream
    auto decl = doc.prepend_child(pugi::node_declaration);
    decl.append_attribute("version")  = "1.0";
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
    for (const auto prefixIter : xsdPrefixesUsed) {
        const auto& xsdType = prefixIter.first;
        const auto& prefix = prefixIter.second;
        if (xsdType == XsdType::NONE || prefix.empty())
            continue;

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
    doc.save(out, "\t", pugi::format_default | pugi::format_no_escapes, pugi::encoding_utf8);
}

void XmlWriter::ToStream(const std::unique_ptr<DataSetBase>& dataset,
                         std::ostream& out)
{ ToStream(*dataset.get(), out); }

} // namespace internal
} // namespace BAM
} // namespace PacBio
