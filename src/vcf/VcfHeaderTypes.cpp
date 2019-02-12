
#include <pbbam/vcf/VcfHeader.h>

#include "pbbam/vcf/VcfHeaderTypes.h"

namespace PacBio {
namespace VCF {

// -------------------
// ContigDefinition
// -------------------

ContigDefinition::ContigDefinition(std::string id) : ContigDefinition(std::move(id), {}) {}

ContigDefinition::ContigDefinition(std::string id,
                                   std::vector<std::pair<std::string, std::string>> attributes)
    : id_{std::move(id)}, attributes_{std::move(attributes)}
{
    if (id_.empty())
        throw std::runtime_error{"VcfFormat: ##contig definition in header has empty ID field"};
}

ContigDefinition::ContigDefinition(const ContigDefinition&) = default;

ContigDefinition::ContigDefinition(ContigDefinition&&) = default;

ContigDefinition& ContigDefinition::operator=(const ContigDefinition&) = default;

ContigDefinition& ContigDefinition::operator=(ContigDefinition&&) = default;

ContigDefinition::~ContigDefinition() = default;

ContigDefinition& ContigDefinition::AddAttribute(std::string id, std::string value)
{
    return AddAttribute(std::make_pair(std::move(id), std::move(value)));
}

ContigDefinition& ContigDefinition::AddAttribute(std::pair<std::string, std::string> attribute)
{
    attributes_.push_back(std::move(attribute));
    return *this;
}

const std::vector<std::pair<std::string, std::string>>& ContigDefinition::Attributes() const
{
    return attributes_;
}

ContigDefinition& ContigDefinition::Attributes(
    std::vector<std::pair<std::string, std::string>> attributes)
{
    attributes_ = std::move(attributes);
    return *this;
}

const std::string& ContigDefinition::Id() const { return id_; }

// -------------------
// FilterDefinition
// -------------------

FilterDefinition::FilterDefinition(std::string id, std::string description)
    : id_{std::move(id)}, description_{std::move(description)}
{
    if (id_.empty())
        throw std::runtime_error{"VcfFormat: FILTER definition in header has empty ID field"};

    if (description_.empty())
        throw std::runtime_error{
            "VcfFormat: FILTER definition in header has empty Description field"};
}

FilterDefinition::FilterDefinition(const FilterDefinition&) = default;

FilterDefinition::FilterDefinition(FilterDefinition&&) = default;

FilterDefinition& FilterDefinition::operator=(const FilterDefinition&) = default;

FilterDefinition& FilterDefinition::operator=(FilterDefinition&&) = default;

FilterDefinition::~FilterDefinition() = default;

const std::string& FilterDefinition::Description() const { return description_; }

const std::string& FilterDefinition::Id() const { return id_; }

// -------------------
// FormatDefinition
// -------------------

FormatDefinition::FormatDefinition(std::string id, std::string number, std::string type,
                                   std::string description)
    : id_{std::move(id)}
    , number_{std::move(number)}
    , type_{std::move(type)}
    , description_{std::move(description)}
{
    if (id_.empty())
        throw std::runtime_error{"VcfFormat: FORMAT definition in header has empty ID field"};

    if (number_.empty())
        throw std::runtime_error{"VcfFormat: FORMAT definition in header has empty Number field"};

    if (type_.empty())
        throw std::runtime_error{"VcfFormat: FORMAT definition in header has empty Type field"};

    if (description_.empty())
        throw std::runtime_error{
            "VcfFormat: FORMAT definition in header has empty Description field"};
}

FormatDefinition::FormatDefinition(const FormatDefinition&) = default;

FormatDefinition::FormatDefinition(FormatDefinition&&) = default;

FormatDefinition& FormatDefinition::operator=(const FormatDefinition&) = default;

FormatDefinition& FormatDefinition::operator=(FormatDefinition&&) = default;

FormatDefinition::~FormatDefinition() = default;

const std::string& FormatDefinition::Description() const { return description_; }

const std::string& FormatDefinition::Id() const { return id_; }

const std::string& FormatDefinition::Number() const { return number_; }

const std::string& FormatDefinition::Type() const { return type_; }

// -------------------
// GeneralDefinition
// -------------------

GeneralDefinition::GeneralDefinition(std::string id, std::string text)
    : id_{std::move(id)}, text_{std::move(text)}
{
    if (id_.empty())
        throw std::runtime_error{
            "VcfFormat: general metadata definition in header has empty label"};

    if (text_.empty())
        throw std::runtime_error{
            "VcfFormat: general metadata definition in header has empty value"};
}

GeneralDefinition::GeneralDefinition(const GeneralDefinition&) = default;

GeneralDefinition::GeneralDefinition(GeneralDefinition&&) = default;

GeneralDefinition& GeneralDefinition::operator=(const GeneralDefinition&) = default;

GeneralDefinition& GeneralDefinition::operator=(GeneralDefinition&&) = default;

GeneralDefinition::~GeneralDefinition() = default;

const std::string& GeneralDefinition::Id() const { return id_; }

const std::string& GeneralDefinition::Text() const { return text_; }

// -------------------
// InfoDefinition
// -------------------

InfoDefinition::InfoDefinition(std::string id, std::string number, std::string type,
                               std::string description, std::string source, std::string version)
    : id_{std::move(id)}
    , number_{std::move(number)}
    , type_{std::move(type)}
    , description_{std::move(description)}
{
    // verify required fields
    if (id_.empty())
        throw std::runtime_error{"VcfFormat: INFO definition in header has empty ID field"};

    if (number_.empty())
        throw std::runtime_error{"VcfFormat: INFO definition in header has empty Number field"};

    if (type_.empty())
        throw std::runtime_error{"VcfFormat: INFO definition in header has empty Type field"};

    if (description_.empty())
        throw std::runtime_error{
            "VcfFormat: INFO definition in header has empty Description field"};

    if (!source.empty()) source_ = std::move(source);
    if (!version.empty()) version_ = std::move(version);
}

InfoDefinition::InfoDefinition(const InfoDefinition&) = default;

InfoDefinition::InfoDefinition(InfoDefinition&&) = default;

InfoDefinition& InfoDefinition::operator=(const InfoDefinition&) = default;

InfoDefinition& InfoDefinition::operator=(InfoDefinition&&) = default;

InfoDefinition::~InfoDefinition() = default;

const std::string& InfoDefinition::Description() const { return description_; }

const std::string& InfoDefinition::Id() const { return id_; }

const std::string& InfoDefinition::Number() const { return number_; }

const boost::optional<std::string>& InfoDefinition::Source() const { return source_; }

InfoDefinition& InfoDefinition::Source(std::string s)
{
    source_ = std::move(s);
    return *this;
}

const std::string& InfoDefinition::Type() const { return type_; }

const boost::optional<std::string>& InfoDefinition::Version() const { return version_; }

InfoDefinition& InfoDefinition::Version(std::string v)
{
    version_ = std::move(v);
    return *this;
}

}  // namespace VCF
}  // namespace PacBio
