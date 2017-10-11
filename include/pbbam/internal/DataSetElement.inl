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

#include "pbbam/internal/DataSetElement.h"

#include <iostream>

namespace PacBio {
namespace BAM {
namespace internal {

// ----------------
// DataSetElement
// ----------------

inline DataSetElement::DataSetElement(const std::string& label, const XsdType& xsd)
    : xsd_(xsd)
    , label_(label)
{ }

inline DataSetElement::DataSetElement(const std::string& label,
                                      const FromInputXml&,
                                      const XsdType& xsd)
    : xsd_(xsd)
    , label_(label, true)
{ }

inline bool DataSetElement::operator==(const DataSetElement& other) const
{
    return xsd_   == other.xsd_   &&
           label_ == other.label_ &&
           text_  == other.text_  &&
           attributes_ == other.attributes_ &&
           children_   == other.children_;
}

inline bool DataSetElement::operator!=(const DataSetElement& other) const
{ return !(*this == other); }

template<typename T>
const T& DataSetElement::operator[](size_t index) const
{ return Child<T>(index); }

template<typename T>
T& DataSetElement::operator[](size_t index)
{ return Child<T>(index); }

template<typename T>
const T& DataSetElement::operator[](const std::string& label) const
{ return Child<T>(label); }

template<typename T>
T& DataSetElement::operator[](const std::string& label)
{ return Child<T>(label); }

inline void DataSetElement::AddChild(const DataSetElement& e)
{ children_.push_back(e); }

inline std::string& DataSetElement::Attribute(const std::string& name)
{ return attributes_[name]; }

inline const std::string& DataSetElement::Attribute(const std::string& name) const
{
    auto iter = attributes_.find(name);
    if (iter == attributes_.cend())
        return SharedNullString();
    return iter->second;
}

inline void DataSetElement::Attribute(const std::string& name, const std::string& value)
{ attributes_[name] = value; }

inline const std::map<std::string, std::string>& DataSetElement::Attributes() const
{ return attributes_; }

inline std::map<std::string, std::string>& DataSetElement::Attributes()
{ return attributes_; }

template<typename T>
inline const T& DataSetElement::Child(size_t index) const
{ return static_cast<const T&>(children_.at(index)); }

template<typename T>
inline T& DataSetElement::Child(size_t index)
{ return static_cast<T&>(children_.at(index)); }

template<typename T>
inline const T& DataSetElement::Child(const std::string& label) const
{ return Child<T>(IndexOf(label)); }

template<typename T>
inline T& DataSetElement::Child(const std::string& label)
{
    const int i = IndexOf(label);
    if (i >= 0) {
        assert(static_cast<size_t>(i) < NumChildren());
        return Child<T>(i);
    } else {
        AddChild(DataSetElement(label));
        return Child<T>(NumChildren()-1);
    }
}

inline const std::vector<DataSetElement>& DataSetElement::Children() const
{ return children_; }

inline std::vector<DataSetElement>& DataSetElement::Children()
{ return children_; }

inline const std::string& DataSetElement::ChildText(const std::string& label) const
{
    if (!HasChild(label))
        return SharedNullString();
    return Child<DataSetElement>(label).Text();
}

inline std::string& DataSetElement::ChildText(const std::string& label)
{
    if (!HasChild(label))
        AddChild(DataSetElement(label));
    return Child<DataSetElement>(label).Text();
}

inline bool DataSetElement::HasAttribute(const std::string& name) const
{ return attributes_.find(name) != attributes_.cend(); }

inline bool DataSetElement::HasChild(const std::string& label) const
{ return IndexOf(label) != -1; }

inline int DataSetElement::IndexOf(const std::string& label) const
{
    const size_t count = NumChildren();
    for (size_t i = 0; i < count; ++i) {
        const DataSetElement& child = children_.at(i);
        if (child.LocalNameLabel() == label || child.label_ == label)
            return i;
    }
    return -1;
}

inline const boost::string_ref DataSetElement::LocalNameLabel() const
{ return label_.LocalName(); }

inline const boost::string_ref DataSetElement::PrefixLabel() const
{ return label_.Prefix(); }

inline const std::string& DataSetElement::QualifiedNameLabel() const
{ return label_.QualifiedName(); }

//inline std::string& DataSetElement::Label()
//{ return label_.QualifiedName(); }

inline void DataSetElement::Label(const std::string& label)
{ label_ = XmlName(label, true); }

inline size_t DataSetElement::NumAttributes() const
{ return attributes_.size(); }

inline size_t DataSetElement::NumChildren() const
{ return children_.size(); }

inline void DataSetElement::RemoveChild(const DataSetElement& e)
{
    children_.erase(
        std::remove(children_.begin(),
                    children_.end(),
                    e),
        children_.end()
    );
}

inline void DataSetElement::ChildText(const std::string& label,
                                         const std::string& text)
{
    if (!HasChild(label)) {
        DataSetElement e(label);
        e.Text(text);
        AddChild(e);
    } else {
        Child<DataSetElement>(label).Text(text);
    }
}

inline bool DataSetElement::IsVerbatimLabel() const
{ return label_.Verbatim(); }

inline const std::string& DataSetElement::Text() const
{ return text_; }

inline std::string& DataSetElement::Text()
{ return text_; }

inline void DataSetElement::Text(const std::string& text)
{ text_ = text; }

inline const XsdType& DataSetElement::Xsd() const
{ return xsd_; }

// ----------------
// XmlName
// ----------------

inline XmlName::XmlName(std::string fullName, bool verbatim)
    : qualifiedName_(std::move(fullName))
    , prefixSize_(0)
    , localNameOffset_(0)
    , localNameSize_(0)
    , verbatim_(verbatim)
{
    const size_t colonFound = qualifiedName_.find(':');
    if (colonFound == std::string::npos || colonFound == 0)
        localNameSize_ = qualifiedName_.size();
    else {
        prefixSize_ = colonFound;
        localNameSize_ = (qualifiedName_.size() - colonFound) - 1;
    }

    // adjust for colon if prefix present
    localNameOffset_ = prefixSize_;
    if (prefixSize_ != 0)
        ++localNameOffset_;
}

inline XmlName::XmlName(const std::string& localName,
                        const std::string& prefix)
    : prefixSize_(prefix.size())
    , localNameOffset_(prefixSize_)
    , localNameSize_(localName.size())
    , verbatim_(true)
{
    qualifiedName_.clear();
    qualifiedName_.reserve(localNameSize_+ prefixSize_ + 1);
    qualifiedName_.append(prefix);
    if (!qualifiedName_.empty())
        qualifiedName_.append(1, ':');
    qualifiedName_.append(localName);

    // adjust for colon if prefix present
    if (prefixSize_ != 0)
        ++localNameOffset_;
}

inline bool XmlName::operator==(const XmlName& other) const
{ return qualifiedName_ == other.qualifiedName_; }

inline bool XmlName::operator!=(const XmlName& other) const
{ return !(*this == other); }

inline const boost::string_ref XmlName::LocalName() const
{ return boost::string_ref(qualifiedName_.data() + localNameOffset_, localNameSize_); }

inline const boost::string_ref XmlName::Prefix() const
{ return boost::string_ref(qualifiedName_.data(), prefixSize_); }

inline const std::string& XmlName::QualifiedName() const
{ return qualifiedName_; }

inline bool XmlName::Verbatim() const
{ return verbatim_; }

} // namespace internal
} // namespace BAM
} // namespace PacBio
