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

#include <algorithm>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
#include <cassert>

namespace PacBio {
namespace BAM {
namespace internal {

class DataSetElement
{
public:
    // allowed to be instantiated directly
    // (for compound objects that don't need a dedicated subtype for children)
    DataSetElement(const std::string& label = std::string());
    DataSetElement(const std::string& label, const std::vector<std::string>& initialChildLabels);
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

    const std::string& Label(void) const;
    std::string& Label(void);

    const std::string& Text(void) const;
    std::string& Text(void);

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
    std::string label_;
    std::string text_;
    std::map<std::string, std::string> attributes_;
    std::vector<DataSetElement> children_;

private:
    int IndexOf(const std::string& label) const;
};

inline DataSetElement::DataSetElement(const std::string& label)
    : label_(label)
{ }

inline DataSetElement::DataSetElement(const std::string& label,
                                      const std::vector<std::string>& initialChildLabels)
    : label_(label)
{
    for (auto childLabel : initialChildLabels)
        AddChild(DataSetElement(childLabel));
}

inline DataSetElement::DataSetElement(const DataSetElement& other)
    : label_(other.label_)
    , text_(other.text_)
    , attributes_(other.attributes_)
    , children_(other.children_)
{ }

inline DataSetElement::DataSetElement(DataSetElement&& other)
    : label_(std::move(other.label_))
    , text_(std::move(other.text_))
    , attributes_(std::move(other.attributes_))
    , children_(std::move(other.children_))
{ }

inline DataSetElement& DataSetElement::operator=(const DataSetElement& other)
{
    label_ = other.label_;
    text_ = other.text_;
    attributes_ = other.attributes_;
    children_ = other.children_;
    return *this;
}

inline DataSetElement& DataSetElement::operator=(DataSetElement&& other)
{
    label_ = std::move(other.label_);
    text_ = std::move(other.text_);
    attributes_ = std::move(other.attributes_);
    children_ = std::move(other.children_);
    return *this;
}

inline DataSetElement::~DataSetElement(void) { }

inline bool DataSetElement::operator==(const DataSetElement& other) const
{
    return label_ == other.label_ &&
           text_  == other.text_ &&
           attributes_ == other.attributes_ &&
           children_ == other.children_;
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

inline const std::map<std::string, std::string>& DataSetElement::Attributes(void) const
{ return attributes_; }

inline std::map<std::string, std::string>& DataSetElement::Attributes(void)
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
        assert(i < NumChildren());
        return Child<T>(i);
    } else {
        AddChild(DataSetElement(label));
        return Child<T>(NumChildren()-1);
    }
}

inline const std::vector<DataSetElement>& DataSetElement::Children(void) const
{ return children_; }

inline std::vector<DataSetElement>& DataSetElement::Children(void)
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
        if (child.label_ == label)
            return i;
    }
    return -1;
}

inline const std::string& DataSetElement::Label(void) const
{ return label_; }

inline std::string& DataSetElement::Label(void)
{ return label_; }

inline void DataSetElement::Label(const std::string& label)
{ label_ = label; }

inline size_t DataSetElement::NumAttributes(void) const
{ return attributes_.size(); }

inline size_t DataSetElement::NumChildren(void) const
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

inline const std::string& DataSetElement::Text(void) const
{ return text_; }

inline std::string& DataSetElement::Text(void)
{ return text_; }

inline void DataSetElement::Text(const std::string& text)
{ text_ = text; }

} // namespace internal
} // namespace BAM
} // namespace PacBio

#endif // DATASETELEMENT_H
