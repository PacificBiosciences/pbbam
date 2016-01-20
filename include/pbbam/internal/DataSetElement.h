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

#ifndef DATASETELEMENT_H
#define DATASETELEMENT_H

#include "pbbam/DataSetXsd.h"

#include <boost/utility/string_ref.hpp>
#include <algorithm>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
#include <cassert>

namespace PacBio {
namespace BAM {
namespace internal {

class XmlName
{
    //    qualified name
    //       |
    //  --------------
    // <pbns:node_name >
    //  ---- ---------
    //   |        |
    //  prefix    local name

public:
    XmlName(const std::string& fullName, bool verbatim = false);
    XmlName(const std::string& localName, const std::string& prefix);
    XmlName(const XmlName& other);
    XmlName(XmlName&& other);
    XmlName& operator=(const XmlName& other);
    XmlName& operator=(XmlName&& other);
    ~XmlName(void);

public:
    bool operator==(const XmlName& other) const;
    bool operator!=(const XmlName& other) const;

public:
    const boost::string_ref LocalName(void) const;
    const boost::string_ref Prefix(void) const;
    const std::string& QualifiedName(void) const;
    bool Verbatim(void) const;

private:
    std::string qualifiedName_;
    size_t prefixSize_;
    size_t localNameOffset_;
    size_t localNameSize_;
    bool verbatim_;
};

struct FromInputXml { };

class DataSetElement
{
public:
    DataSetElement(const std::string& label, const XsdType& xsd = XsdType::NONE);
    DataSetElement(const std::string& label, const FromInputXml& fromInputXml, const XsdType& xsd = XsdType::NONE);
    DataSetElement(const DataSetElement& other);
    DataSetElement(DataSetElement&& other);
    DataSetElement& operator=(const DataSetElement& other);
    DataSetElement& operator=(DataSetElement&& other);
    virtual ~DataSetElement(void);

public:
    bool operator==(const DataSetElement& other) const;
    bool operator!=(const DataSetElement& other) const;

public:
    const std::string& Attribute(const std::string& name) const;
    std::string& Attribute(const std::string& name);
    const std::map<std::string, std::string>& Attributes(void) const;
    std::map<std::string, std::string>& Attributes(void);
    bool HasAttribute(const std::string& name) const;

    const std::vector<DataSetElement>& Children(void) const;
    std::vector<DataSetElement>& Children(void);
    bool HasChild(const std::string& label) const;

    const boost::string_ref LocalNameLabel(void) const;
    const boost::string_ref PrefixLabel(void) const;
    const std::string& QualifiedNameLabel(void) const;
    bool IsVerbatimLabel(void) const;

    const std::string& Text(void) const;
    std::string& Text(void);

    const XsdType& Xsd(void) const;

public:
    void Attribute(const std::string& name, const std::string& value);
    void Label(const std::string& label);
    void Text(const std::string& text);

public:
    size_t NumAttributes(void) const;
    size_t NumChildren(void) const;

public:
    void AddChild(const DataSetElement& e);
    void RemoveChild(const DataSetElement& e);

    template<typename T>
    const T& Child(size_t index) const;

    template<typename T>
    T& Child(size_t index);

    template<typename T>
    const T& Child(const std::string& label) const;

    template<typename T>
    T& Child(const std::string& label);

    template<typename T>
    const T& operator[](size_t index) const;

    template<typename T>
    T& operator[](size_t index);

    template<typename T = DataSetElement>
    const T& operator[](const std::string& label) const;

    template<typename T = DataSetElement>
    T& operator[](const std::string& label);

protected:
    static const std::string& SharedNullString(void);

public:
    const std::string& ChildText(const std::string& label) const;
    std::string& ChildText(const std::string& label);
    void ChildText(const std::string& label, const std::string& text);

protected:
    XsdType xsd_;
    XmlName label_;
    std::string text_;
    std::map<std::string, std::string> attributes_;
    std::vector<DataSetElement> children_;

private:
    int IndexOf(const std::string& label) const;
};

} // namespace internal
} // namespace BAM
} // namespace PacBio

#include "pbbam/internal/DataSetElement.inl"

#endif // DATASETELEMENT_H
