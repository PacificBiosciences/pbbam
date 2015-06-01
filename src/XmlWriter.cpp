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
#include "pbbam/dataset/DataSetBase.h"
#include "pugixml/pugixml.hpp"
#include <fstream>
#include <iostream>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

static
void ToXml(const DataSetElement& node,
           pugi::xml_node& parentXml)
{
    // create child of parent, w/ label & text
    if ( node.Label().empty() )
        return; // error?
    pugi::xml_node xmlNode = parentXml.append_child(node.Label().c_str());

    if (!node.Text().empty())
        xmlNode.text().set(node.Text().c_str());

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
        ToXml(child, xmlNode);
    }
}

} // namespace internal
} // namespace BAM
} // namespace PacBio

void XmlWriter::ToStream(const DataSetBase& dataset, std::ostream& out)
{
    pugi::xml_document doc;

    // create top-level dataset XML node
    const string& label = dataset.Label();
    if (label.empty())
        throw std::exception();
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

    // iterate children, recursively building up subtree
    auto childIter = dataset.Children().cbegin();
    auto childEnd  = dataset.Children().cend();
    for ( ; childIter != childEnd; ++childIter) {
        const DataSetElement& child = (*childIter);
        ToXml(child, root);
    }

    // write XML to stream
    pugi::xml_node decl = doc.prepend_child(pugi::node_declaration);
    decl.append_attribute("version")  = "1.0";
    decl.append_attribute("encoding") = "utf-8";

    // "no escapes" to allow explicit ">" "<" comparison operators in filter parameters
    // we may remove this if/when comparison is separated from the value
    doc.save(out, "\t", pugi::format_default | pugi::format_no_escapes, pugi::encoding_utf8);
}
