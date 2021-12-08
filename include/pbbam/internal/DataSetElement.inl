#ifndef PBBAM_DATASETELEMENT_INL
#define PBBAM_DATASETELEMENT_INL

#include <pbbam/Config.h>

#include <pbbam/internal/DataSetElement.h>

#include <iostream>
#include <stdexcept>
#include <tuple>
#include <typeinfo>

#include <cassert>

namespace PacBio {
namespace BAM {
namespace internal {

// ----------------
// DataSetElement
// ----------------

inline DataSetElement::DataSetElement(const std::string& label, const XsdType& xsd)
    : xsd_(xsd), label_(label)
{
}

inline DataSetElement::DataSetElement(const std::string& label, const FromInputXml&,
                                      const XsdType& xsd)
    : xsd_(xsd), label_(label, true)
{
}

inline bool DataSetElement::operator==(const DataSetElement& other) const noexcept
{
    return std::tie(xsd_, label_, text_, attributes_, children_) ==
           std::tie(other.xsd_, other.label_, other.text_, other.attributes_, other.children_);
}

inline bool DataSetElement::operator!=(const DataSetElement& other) const noexcept
{
    return !(*this == other);
}

template <typename T>
const T& DataSetElement::operator[](size_t index) const
{
    return Child<T>(index);
}

template <typename T>
T& DataSetElement::operator[](size_t index)
{
    return Child<T>(index);
}

template <typename T>
const T& DataSetElement::operator[](const std::string& label) const
{
    return Child<T>(label);
}

template <typename T>
T& DataSetElement::operator[](const std::string& label)
{
    return Child<T>(label);
}

template <typename T>
void DataSetElement::AddChild(const T& e)
{
    children_.push_back(std::make_shared<T>(e));
}

inline void DataSetElement::AddChild(const DataSetElement& e)
{
    children_.push_back(std::make_shared<DataSetElement>(e));
}

inline void DataSetElement::AddChild(std::shared_ptr<DataSetElement> e) { children_.push_back(e); }

inline std::string& DataSetElement::Attribute(const std::string& name) { return attributes_[name]; }

inline const std::string& DataSetElement::Attribute(const std::string& name) const
{
    auto iter = attributes_.find(name);
    if (iter == attributes_.cend()) {
        return SharedNullString();
    }
    return iter->second;
}

inline void DataSetElement::Attribute(const std::string& name, const std::string& value)
{
    attributes_[name] = value;
}

inline const std::map<std::string, std::string>& DataSetElement::Attributes() const
{
    return attributes_;
}

inline std::map<std::string, std::string>& DataSetElement::Attributes() { return attributes_; }

template <typename T>
const T& DataSetElement::Child(size_t index) const
{
    DataSetElement* child = children_.at(index).get();
    if (child == nullptr) {
        throw std::runtime_error{
            "[pbbam] dataset element ERROR: cannot access null child at index " +
            std::to_string(index) + " in element: " + QualifiedNameLabel()};
    }
    const T* c = dynamic_cast<const T*>(child);
    return *c;
}

template <typename T>
T& DataSetElement::Child(size_t index)
{
    DataSetElement* child = children_.at(index).get();
    if (child == nullptr) {
        throw std::runtime_error{
            "[pbbam] dataset element ERROR: cannot access null child at index " +
            std::to_string(index) + " in element: " + QualifiedNameLabel()};
    }
    T* c = dynamic_cast<T*>(child);
    return *c;
}

template <typename T>
const T& DataSetElement::Child(const std::string& label) const
{
    const auto index = IndexOf(label);
    return Child<T>(index);
}

template <typename T>
T& DataSetElement::Child(const std::string& label)
{
    const int i = IndexOf(label);
    if (i >= 0) {
        assert(static_cast<size_t>(i) < NumChildren());
        return Child<T>(i);
    } else {
        AddChild(T());
        return Child<T>(NumChildren() - 1);
    }
}

template <>
inline DataSetElement& DataSetElement::Child<DataSetElement>(const std::string& label)
{
    const int i = IndexOf(label);
    if (i >= 0) {
        assert(static_cast<size_t>(i) < NumChildren());
        return Child<DataSetElement>(i);
    } else {
        AddChild(DataSetElement{label});
        return Child<DataSetElement>(NumChildren() - 1);
    }
}

inline const std::vector<std::shared_ptr<DataSetElement>>& DataSetElement::Children() const
{
    return children_;
}

inline std::vector<std::shared_ptr<DataSetElement>>& DataSetElement::Children()
{
    return children_;
}

inline const std::string& DataSetElement::ChildText(const std::string& label) const
{
    if (!HasChild(label)) {
        return SharedNullString();
    }
    return Child<DataSetElement>(label).Text();
}

inline std::string& DataSetElement::ChildText(const std::string& label)
{
    if (!HasChild(label)) {
        AddChild(DataSetElement(label));
    }
    return Child<DataSetElement>(label).Text();
}

inline bool DataSetElement::HasAttribute(const std::string& name) const
{
    return attributes_.find(name) != attributes_.cend();
}

inline bool DataSetElement::HasChild(const std::string& label) const
{
    return IndexOf(label) != -1;
}

inline int DataSetElement::IndexOf(const std::string& label) const
{
    const size_t count = NumChildren();
    for (size_t i = 0; i < count; ++i) {
        const DataSetElement& child = *(children_.at(i).get());
        if (child.LocalNameLabel() == label || child.QualifiedNameLabel() == label ||
            child.label_ == label) {
            return i;
        }
    }
    return -1;
}

inline const boost::string_ref DataSetElement::LocalNameLabel() const { return label_.LocalName(); }

inline const boost::string_ref DataSetElement::PrefixLabel() const { return label_.Prefix(); }

inline const std::string& DataSetElement::QualifiedNameLabel() const
{
    return label_.QualifiedName();
}

inline void DataSetElement::Label(const std::string& label) { label_ = XmlName(label, true); }

inline size_t DataSetElement::NumAttributes() const { return attributes_.size(); }

inline size_t DataSetElement::NumChildren() const { return children_.size(); }

inline size_t DataSetElement::Size() const { return children_.size(); }

inline void DataSetElement::RemoveChild(const DataSetElement& e)
{
    std::vector<std::shared_ptr<DataSetElement>> newChildren;
    for (std::shared_ptr<DataSetElement>& child : children_) {
        if (*(child.get()) != e) {
            newChildren.push_back(std::move(child));
        }
    }
    children_ = std::move(newChildren);
}

inline void DataSetElement::ChildText(const std::string& label, const std::string& text)
{
    if (!HasChild(label)) {
        DataSetElement e(label);
        e.Text(text);
        AddChild(e);
    } else {
        Child<DataSetElement>(label).Text(text);
    }
}

inline bool DataSetElement::IsVerbatimLabel() const { return label_.Verbatim(); }

inline const std::string& DataSetElement::Text() const { return text_; }

inline std::string& DataSetElement::Text() { return text_; }

inline void DataSetElement::Text(const std::string& text) { text_ = text; }

inline const XsdType& DataSetElement::Xsd() const { return xsd_; }

// ----------------------------
// DataSetElementIteratorBase
// ----------------------------

inline DataSetElementIteratorBase::DataSetElementIteratorBase(const DataSetElement* parent,
                                                              size_t i)
    : parent_(parent), index_(i)
{
}

inline bool DataSetElementIteratorBase::operator==(const DataSetElementIteratorBase& other) const
    noexcept
{
    return std::tie(parent_, index_) == std::tie(other.parent_, other.index_);
}

inline bool DataSetElementIteratorBase::operator!=(const DataSetElementIteratorBase& other) const
    noexcept
{
    return !(*this == other);
}

inline void DataSetElementIteratorBase::Advance()
{
    if (index_ >= parent_->NumChildren()) {
        parent_ = nullptr;
        return;
    }
    ++index_;
}

// ------------------------
// DataSetElementIterator
// ------------------------

template <typename T>
DataSetElementIterator<T>::DataSetElementIterator(const DataSetElement* parent, size_t i)
    : DataSetElementIteratorBase(parent, i)
{
}

template <typename T>
T& DataSetElementIterator<T>::operator*() noexcept
{
    return parent_->template Child<T>(index_);
}

template <typename T>
T* DataSetElementIterator<T>::operator->() noexcept
{
    return &(operator*());
}

template <typename T>
DataSetElementIterator<T>& DataSetElementIterator<T>::operator++()
{
    Advance();
    return *this;
}

template <typename T>
DataSetElementIterator<T> DataSetElementIterator<T>::operator++(int)
{
    DataSetElementIterator<T> result(*this);
    ++(*this);
    return result;
}

// -----------------------------
// DataSetElementConstIterator
// -----------------------------

template <typename T>
DataSetElementConstIterator<T>::DataSetElementConstIterator(const DataSetElement* parent, size_t i)
    : DataSetElementIteratorBase(parent, i)
{
}

template <typename T>
const T& DataSetElementConstIterator<T>::operator*() const noexcept
{
    return parent_->template Child<const T>(index_);
}

template <typename T>
const T* DataSetElementConstIterator<T>::operator->() const noexcept
{
    return &(operator*());
}

template <typename T>
DataSetElementConstIterator<T>& DataSetElementConstIterator<T>::operator++()
{
    Advance();
    return *this;
}

template <typename T>
DataSetElementConstIterator<T> DataSetElementConstIterator<T>::operator++(int)
{
    DataSetElementConstIterator<T> result(*this);
    ++(*this);
    return result;
}

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
    if (colonFound == std::string::npos || colonFound == 0) {
        localNameSize_ = qualifiedName_.size();
    } else {
        prefixSize_ = colonFound;
        localNameSize_ = (qualifiedName_.size() - colonFound) - 1;
    }

    // adjust for colon if prefix present
    localNameOffset_ = prefixSize_;
    if (prefixSize_ != 0) {
        ++localNameOffset_;
    }
}

inline XmlName::XmlName(const std::string& localName, const std::string& prefix)
    : prefixSize_(prefix.size())
    , localNameOffset_(prefixSize_)
    , localNameSize_(localName.size())
    , verbatim_(true)
{
    qualifiedName_.clear();
    qualifiedName_.reserve(localNameSize_ + prefixSize_ + 1);
    qualifiedName_.append(prefix);
    if (!qualifiedName_.empty()) {
        qualifiedName_.append(1, ':');
    }
    qualifiedName_.append(localName);

    // adjust for colon if prefix present
    if (prefixSize_ != 0) {
        ++localNameOffset_;
    }
}

inline bool XmlName::operator==(const XmlName& other) const noexcept
{
    return qualifiedName_ == other.qualifiedName_;
}

inline bool XmlName::operator!=(const XmlName& other) const noexcept { return !(*this == other); }

inline const boost::string_ref XmlName::LocalName() const
{
    return boost::string_ref(qualifiedName_.data() + localNameOffset_, localNameSize_);
}

inline const boost::string_ref XmlName::Prefix() const
{
    return boost::string_ref(qualifiedName_.data(), prefixSize_);
}

inline const std::string& XmlName::QualifiedName() const { return qualifiedName_; }

inline bool XmlName::Verbatim() const { return verbatim_; }

}  // namespace internal
}  // namespace BAM
}  // namespace PacBio

#endif // PBBAM_DATASETELEMENT_INL
