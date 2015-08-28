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

#include "XmlReader.h"
#include "StringUtils.h"
#include "pugixml/pugixml.hpp"
#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include <cassert>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

static
void UpdateRegistry(const string& attributeName,
                    const string& attributeValue,
                    NamespaceRegistry& registry)
{
    vector<string> nameParts = Split(attributeName, ':');
    assert(!nameParts.empty());
    if (nameParts.size() > 2)
        throw std::runtime_error("malformed xmlns attribute: " + attributeName);

    const bool isDefault = (nameParts.size() == 1);
    const XsdType& xsd = registry.XsdForUri(attributeValue);

    if (isDefault)
        registry.SetDefaultXsd(xsd);
    else {
        assert(nameParts.size() == 2);
        const string& name = nameParts.at(1);
        const string& uri  = attributeValue;
        NamespaceInfo namespaceInfo(name, uri);
        registry.Register(xsd, namespaceInfo);
    }
}

static
void FromXml(const pugi::xml_node& xmlNode, DataSetElement& parent)
{
    // ignore non-named XML nodes
    //
    // pugi::xml separates XML parts into more node types than we use
    //
    const string& label = xmlNode.name();
    if (label.empty())
        return;

    // label & text
    DataSetElement e(xmlNode.name(), FromInputXml());
    e.Text(xmlNode.text().get());

    // iterate attributes
    auto attrIter = xmlNode.attributes_begin();
    auto attrEnd  = xmlNode.attributes_end();
    for ( ; attrIter != attrEnd; ++attrIter )
        e.Attribute(attrIter->name(), attrIter->value());

    // iterate children, recursively building up subtree
    auto childIter = xmlNode.begin();
    auto childEnd = xmlNode.end();
    for ( ; childIter != childEnd; ++childIter ) {
        pugi::xml_node childNode = *childIter;
        FromXml(childNode, e);
    }

    // add our element to its parent
    parent.AddChild(e);
}

} // namespace internal
} // namespace BAM
} // namespace PacBio

std::unique_ptr<DataSetBase> XmlReader::FromStream(istream& in)
{
    pugi::xml_document doc;
    const pugi::xml_parse_result& loadResult = doc.load(in);
    if (loadResult.status != pugi::status_ok)
        throw std::runtime_error(string("could not read XML file, error code:") + to_string(loadResult.status) );

    // parse top-level attributes
    pugi::xml_node rootNode = doc.document_element();
    if (rootNode == pugi::xml_node())
        throw std::runtime_error("could not fetch XML root node");

    // create dataset matching type strings
    std::unique_ptr<DataSetBase> dataset(new DataSetBase);
    dataset->Label(rootNode.name());

    // iterate attributes, capture namespace info
    const string xmlnsPrefix("xmlns");
    auto attrIter = rootNode.attributes_begin();
    auto attrEnd  = rootNode.attributes_end();
    for ( ; attrIter != attrEnd; ++attrIter ) {
        const string& name = attrIter->name();
        const string& value = attrIter->value();
        dataset->Attribute(name, value);

        if (name.find(xmlnsPrefix) == 0)
            UpdateRegistry(name, value, dataset->Namespaces());
    }

    // iterate children, recursively building up subtree
    auto childIter = rootNode.begin();
    auto childEnd = rootNode.end();
    for ( ; childIter != childEnd; ++childIter ) {
        pugi::xml_node childNode = *childIter;
        internal::FromXml(childNode, *dataset.get());
    }

    return dataset;
}
