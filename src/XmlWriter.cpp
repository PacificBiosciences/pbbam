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

#include "XmlWriter.h"
#include "pbbam/DataSet.h"
#include "pugixml/pugixml.hpp"
#include <fstream>
#include <iostream>
#include <map>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

static
string Prefix(const string& input)
{
    const size_t colonFound = input.find(':');
    if (colonFound == std::string::npos || colonFound == 0)
        return string();
    return input.substr(0, colonFound);
}

static
string OutputName(const DataSetElement& node,
                  const NamespaceRegistry& registry)
{
    // if from input XML, respect the namespaces given
    if (node.IsVerbatimLabel()) 
        return node.QualifiedNameLabel();

    // otherwise, probably user-generated
    else {
        // if no namespace prefix, prepend the appropriate one & return
        if (node.PrefixLabel().empty()) {
            static const string colon = ":";
            XsdType xsdType = node.Xsd();
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
           map<XsdType, string>& xsdPrefixesUsed,
           pugi::xml_node& parentXml)
{
    // create child of parent, w/ label & text
    const string& label = OutputName(node, registry);
    if (label.empty())
        return; // error?
    pugi::xml_node xmlNode = parentXml.append_child(label.c_str());

    if (!node.Text().empty())
        xmlNode.text().set(node.Text().c_str());

    // store XSD type for later
    const string prefix = Prefix(label);
    if (!prefix.empty())
        xsdPrefixesUsed[node.Xsd()] = prefix;

    // add attributes
    auto attrIter = node.Attributes().cbegin();
    auto attrEnd  = node.Attributes().cend();
    for ( ; attrIter != attrEnd; ++attrIter) {
        const string& name = attrIter->first;
        if (name.empty())
            continue;
        pugi::xml_attribute attr = xmlNode.append_attribute(name.c_str());
        attr.set_value(attrIter->second.c_str());
    }

    // additional stuff later? (e.g. comments)

    // iterate children, recursively building up subtree
    auto childIter = node.Children().cbegin();
    auto childEnd  = node.Children().cend();
    for ( ; childIter != childEnd; ++childIter) {
        const DataSetElement& child = (*childIter);
        ToXml(child, registry, xsdPrefixesUsed, xmlNode);
    }
}

} // namespace internal
} // namespace BAM
} // namespace PacBio

void XmlWriter::ToStream(const DataSetBase& dataset,
                         ostream& out)
{
    pugi::xml_document doc;

    const NamespaceRegistry& registry = dataset.Namespaces();

    // create top-level dataset XML node
    const string& label = internal::OutputName(dataset, registry);
    if (label.empty())
        throw std::runtime_error("could not convert dataset node to XML");
    pugi::xml_node root = doc.append_child(label.c_str());

    const string& text = dataset.Text();
    if (!text.empty())
        root.text().set(text.c_str());

    // add top-level attributes
    auto attrIter = dataset.Attributes().cbegin();
    auto attrEnd  = dataset.Attributes().cend();
    for ( ; attrIter != attrEnd; ++attrIter) {
        const string name = attrIter->first;
        const string value = attrIter->second;
        if (name.empty())
            continue;
        pugi::xml_attribute attr = root.append_attribute(name.c_str());
        attr.set_value(value.c_str());
    }

    map<XsdType, string> xsdPrefixesUsed;
    xsdPrefixesUsed[dataset.Xsd()] = Prefix(label);

    // iterate children, recursively building up subtree
    auto childIter = dataset.Children().cbegin();
    auto childEnd  = dataset.Children().cend();
    for ( ; childIter != childEnd; ++childIter) {
        const DataSetElement& child = (*childIter);
        ToXml(child, registry, xsdPrefixesUsed, root);
    }

    // write XML to stream
    pugi::xml_node decl = doc.prepend_child(pugi::node_declaration);
    decl.append_attribute("version")  = "1.0";
    decl.append_attribute("encoding") = "utf-8";

    // add XSD namespace attributes
    pugi::xml_attribute xmlnsDefaultAttribute = root.attribute("xmlns");
    if (xmlnsDefaultAttribute.empty()) {
        xmlnsDefaultAttribute = root.append_attribute("xmlns");
        xmlnsDefaultAttribute.set_value(registry.DefaultNamespace().Uri().c_str());
    }
    pugi::xml_attribute xsiAttribute = root.attribute("xmlns:xsi");
    if (xsiAttribute.empty()) {
        xsiAttribute = root.append_attribute("xmlns:xsi");
        xsiAttribute.set_value("http://www.w3.org/2001/XMLSchema-instance");
    }
    pugi::xml_attribute xsiSchemaLocationAttribute = root.attribute("xsi:schemaLocation");
    if (xsiSchemaLocationAttribute.empty()) {
        xsiSchemaLocationAttribute = root.append_attribute("xsi:schemaLocation");
        xsiSchemaLocationAttribute.set_value(registry.DefaultNamespace().Uri().c_str());
    }

    static const string xmlnsPrefix = "xmlns:";
    map<XsdType, string>::const_iterator prefixIter = xsdPrefixesUsed.cbegin();
    map<XsdType, string>::const_iterator prefixEnd  = xsdPrefixesUsed.cend();
    for ( ; prefixIter != prefixEnd; ++prefixIter ) {
        const XsdType& xsd = prefixIter->first;
        const string& prefix = prefixIter->second;
        if (xsd == XsdType::NONE || prefix.empty())
            continue;
        const NamespaceInfo& nsInfo = registry.Namespace(xsd);
        assert(nsInfo.Name() == prefix);
        const string xmlnsName = xmlnsPrefix + prefix;
        pugi::xml_attribute xmlnsAttribute = root.attribute(xmlnsName.c_str());
        if (xmlnsAttribute.empty()) {
            xmlnsAttribute = root.append_attribute(xmlnsName.c_str());
            xmlnsAttribute.set_value(nsInfo.Uri().c_str());
        }
    }

    // "no escapes" to allow explicit ">" "<" comparison operators in filter parameters
    // we may remove this if/when comparison is separated from the value
    doc.save(out, "\t", pugi::format_default | pugi::format_no_escapes, pugi::encoding_utf8);
}

void XmlWriter::ToStream(const unique_ptr<DataSetBase>& dataset,
                         ostream& out)
{ ToStream(*dataset.get(), out); }
