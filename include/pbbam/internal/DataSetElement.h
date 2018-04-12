// Author: Derek Barnett

#ifndef DATASETELEMENT_H
#define DATASETELEMENT_H

#include "pbbam/DataSetXsd.h"

#include <algorithm>
#include <boost/utility/string_ref.hpp>
#include <cassert>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

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
    XmlName(std::string fullName, bool verbatim = false);
    XmlName(const std::string& localName, const std::string& prefix);

    XmlName(const XmlName&) = default;
    XmlName(XmlName&&) = default;
    XmlName& operator=(const XmlName&) = default;
    XmlName& operator=(XmlName&&) = default;
    ~XmlName() = default;

public:
    bool operator==(const XmlName& other) const;
    bool operator!=(const XmlName& other) const;

public:
    const boost::string_ref LocalName() const;
    const boost::string_ref Prefix() const;
    const std::string& QualifiedName() const;
    bool Verbatim() const;

private:
    std::string qualifiedName_;
    size_t prefixSize_;
    size_t localNameOffset_;
    size_t localNameSize_;
    bool verbatim_;
};

struct FromInputXml
{
};

class DataSetElement
{
public:
    DataSetElement(const std::string& label, const XsdType& xsd = XsdType::NONE);
    DataSetElement(const std::string& label, const FromInputXml& fromInputXml,
                   const XsdType& xsd = XsdType::NONE);

    DataSetElement(const DataSetElement&) = default;
    DataSetElement(DataSetElement&&) = default;
    DataSetElement& operator=(const DataSetElement&) = default;
    DataSetElement& operator=(DataSetElement&&) = default;
    virtual ~DataSetElement() = default;

public:
    bool operator==(const DataSetElement& other) const;
    bool operator!=(const DataSetElement& other) const;

public:
    const std::string& Attribute(const std::string& name) const;
    std::string& Attribute(const std::string& name);
    const std::map<std::string, std::string>& Attributes() const;
    std::map<std::string, std::string>& Attributes();
    bool HasAttribute(const std::string& name) const;

    const std::vector<DataSetElement>& Children() const;
    std::vector<DataSetElement>& Children();
    bool HasChild(const std::string& label) const;

    const boost::string_ref LocalNameLabel() const;
    const boost::string_ref PrefixLabel() const;
    const std::string& QualifiedNameLabel() const;
    bool IsVerbatimLabel() const;

    const std::string& Text() const;
    std::string& Text();

    const XsdType& Xsd() const;

public:
    void Attribute(const std::string& name, const std::string& value);
    void Label(const std::string& label);
    void Text(const std::string& text);

public:
    size_t NumAttributes() const;
    size_t NumChildren() const;

public:
    void AddChild(const DataSetElement& e);
    void RemoveChild(const DataSetElement& e);

    template <typename T>
    const T& Child(size_t index) const;

    template <typename T>
    T& Child(size_t index);

    template <typename T>
    const T& Child(const std::string& label) const;

    template <typename T>
    T& Child(const std::string& label);

    template <typename T>
    const T& operator[](size_t index) const;

    template <typename T>
    T& operator[](size_t index);

    template <typename T = DataSetElement>
    const T& operator[](const std::string& label) const;

    template <typename T = DataSetElement>
    T& operator[](const std::string& label);

protected:
    static const std::string& SharedNullString();

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

}  // namespace internal
}  // namespace BAM
}  // namespace PacBio

#include "pbbam/internal/DataSetElement.inl"

#endif  // DATASETELEMENT_H
