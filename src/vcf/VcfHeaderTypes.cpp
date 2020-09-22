#include "../PbbamInternalConfig.h"

#include <pbbam/vcf/VcfHeaderTypes.h>

#include <cassert>

#include <type_traits>

#include <pbbam/vcf/VcfHeader.h>

#include "VcfFormatException.h"

namespace PacBio {
namespace VCF {

// -------------------
// ContigDefinition
// -------------------

static_assert(std::is_copy_constructible<ContigDefinition>::value,
              "ContigDefinition(const ContigDefinition&) is not = default");
static_assert(std::is_copy_assignable<ContigDefinition>::value,
              "ContigDefinition& operator=(const ContigDefinition&) is not = default");

static_assert(std::is_nothrow_move_constructible<ContigDefinition>::value,
              "ContigDefinition(ContigDefinition&&) is not = noexcept");
static_assert(std::is_nothrow_move_assignable<ContigDefinition>::value ==
                  std::is_nothrow_move_assignable<std::string>::value,
              "");

ContigDefinition::ContigDefinition(std::string id) : ContigDefinition{std::move(id), {}} {}

// clang-format off
ContigDefinition::ContigDefinition(std::string id,
                                   std::vector<std::pair<std::string, std::string>> attributes)
    : id_{std::move(id)}, attributes_{std::move(attributes)}
{
    if (id_.empty())
        throw VcfFormatException{"##contig definition in header has empty ID field"};
}
// clang-format on

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

static_assert(std::is_copy_constructible<FilterDefinition>::value,
              "FilterDefinition(const FilterDefinition&) is not = default");
static_assert(std::is_copy_assignable<FilterDefinition>::value,
              "FilterDefinition& operator=(const FilterDefinition&) is not = default");

static_assert(std::is_nothrow_move_constructible<FilterDefinition>::value,
              "FilterDefinition(FilterDefinition&&) is not = noexcept");
static_assert(std::is_nothrow_move_assignable<FilterDefinition>::value ==
                  std::is_nothrow_move_assignable<std::string>::value,
              "");

// clang-format off
FilterDefinition::FilterDefinition(std::string id, std::string description)
    : id_{std::move(id)}, description_{std::move(description)}
{
    if (id_.empty())
        throw VcfFormatException{"FILTER definition in header has empty ID field"};

    if (description_.empty())
        throw VcfFormatException{"FILTER definition in header has empty Description field"};
}
// clang-format on

const std::string& FilterDefinition::Description() const { return description_; }

const std::string& FilterDefinition::Id() const { return id_; }

// -------------------
// FormatDefinition
// -------------------

static_assert(std::is_copy_constructible<FormatDefinition>::value,
              "FormatDefinition(const FormatDefinition&) is not = default");
static_assert(std::is_copy_assignable<FormatDefinition>::value,
              "FormatDefinition& operator=(const FormatDefinition&) is not = default");

static_assert(std::is_nothrow_move_constructible<FormatDefinition>::value,
              "FormatDefinition(FormatDefinition&&) is not = noexcept");
static_assert(std::is_nothrow_move_assignable<FormatDefinition>::value ==
                  std::is_nothrow_move_assignable<std::string>::value,
              "");

// clang-format off
FormatDefinition::FormatDefinition(std::string id, std::string number, std::string type,
                                   std::string description)
    : id_{std::move(id)}
    , number_{std::move(number)}
    , type_{std::move(type)}
    , description_{std::move(description)}
{
    if (id_.empty())
        throw VcfFormatException{"FORMAT definition in header has empty ID field"};

    if (number_.empty())
        throw VcfFormatException{"FORMAT definition in header has empty Number field"};

    if (type_.empty())
        throw VcfFormatException{"FORMAT definition in header has empty Type field"};

    if (description_.empty())
        throw VcfFormatException{"FORMAT definition in header has empty Description field"};
}
// clang-format on

const std::string& FormatDefinition::Description() const { return description_; }

const std::string& FormatDefinition::Id() const { return id_; }

const std::string& FormatDefinition::Number() const { return number_; }

const std::string& FormatDefinition::Type() const { return type_; }

// -------------------
// GeneralDefinition
// -------------------

static_assert(std::is_copy_constructible<GeneralDefinition>::value,
              "GeneralDefinition(const GeneralDefinition&) is not = default");
static_assert(std::is_copy_assignable<GeneralDefinition>::value,
              "GeneralDefinition& operator=(const GeneralDefinition&) is not = default");

static_assert(std::is_nothrow_move_constructible<GeneralDefinition>::value,
              "GeneralDefinition(GeneralDefinition&&) is not = noexcept");
static_assert(std::is_nothrow_move_assignable<GeneralDefinition>::value ==
                  std::is_nothrow_move_assignable<std::string>::value,
              "");

// clang-format off
GeneralDefinition::GeneralDefinition(std::string id, std::string text)
    : id_{std::move(id)}, text_{std::move(text)}
{
    if (id_.empty())
        throw VcfFormatException{"general metadata definition in header has empty label"};

    if (text_.empty())
        throw VcfFormatException{"general metadata definition in header has empty value"};
}
// clang-format on

const std::string& GeneralDefinition::Id() const { return id_; }

const std::string& GeneralDefinition::Text() const { return text_; }

// -------------------
// InfoDefinition
// -------------------

static_assert(std::is_copy_constructible<InfoDefinition>::value,
              "InfoDefinition(const InfoDefinition&) is not = default");
static_assert(std::is_copy_assignable<InfoDefinition>::value,
              "InfoDefinition& operator=(const InfoDefinition&) is not = default");

#ifndef __INTEL_COMPILER
static_assert(std::is_nothrow_move_constructible<InfoDefinition>::value,
              "InfoDefinition(InfoDefinition&&) is not = noexcept");
static_assert(std::is_nothrow_move_assignable<InfoDefinition>::value ==
                  std::is_nothrow_move_assignable<std::string>::value,
              "");
#endif

// clang-format off
InfoDefinition::InfoDefinition(std::string id, std::string number, std::string type,
                               std::string description, std::string source, std::string version)
    : id_{std::move(id)}
    , number_{std::move(number)}
    , type_{std::move(type)}
    , description_{std::move(description)}
{
    // verify required fields
    if (id_.empty())
        throw VcfFormatException{"INFO definition in header has empty ID field"};

    if (number_.empty())
        throw VcfFormatException{"INFO definition in header has empty Number field"};

    if (type_.empty())
        throw VcfFormatException{"INFO definition in header has empty Type field"};

    if (description_.empty())
        throw VcfFormatException{"INFO definition in header has empty Description field"};

    if (!source.empty()) source_ = std::move(source);
    if (!version.empty()) version_ = std::move(version);
}
// clang-format on

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
