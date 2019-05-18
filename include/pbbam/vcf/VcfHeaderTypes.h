
#ifndef PBBAM_VCF_VCFHEADERTYPES_H
#define PBBAM_VCF_VCFHEADERTYPES_H

#include <string>
#include <utility>
#include <vector>

#include <boost/optional.hpp>

#include "pbbam/Config.h"

namespace PacBio {
namespace VCF {

using Sample = std::string;

class ContigDefinition
{
public:
    explicit ContigDefinition(std::string id);
    ContigDefinition(std::string id, std::vector<std::pair<std::string, std::string>> attributes);

    ContigDefinition() = delete;
    ContigDefinition(const ContigDefinition&);
    ContigDefinition(ContigDefinition&&) noexcept;
    ContigDefinition& operator=(const ContigDefinition&);
    ContigDefinition& operator=(ContigDefinition&&) PBBAM_NOEXCEPT_MOVE_ASSIGN;
    ~ContigDefinition();

public:
    const std::string& Id() const;
    const std::vector<std::pair<std::string, std::string>>& Attributes() const;

    ContigDefinition& AddAttribute(std::string id, std::string value);
    ContigDefinition& AddAttribute(std::pair<std::string, std::string> attribute);
    ContigDefinition& Attributes(std::vector<std::pair<std::string, std::string>> attributes);

private:
    std::string id_;
    std::vector<std::pair<std::string, std::string>> attributes_;
};

///
/// \brief The FilterDefinition class
///
class FilterDefinition
{
public:
    FilterDefinition(std::string id, std::string description);

    FilterDefinition() = delete;
    FilterDefinition(const FilterDefinition&);
    FilterDefinition(FilterDefinition&&) noexcept;
    FilterDefinition& operator=(const FilterDefinition&);
    FilterDefinition& operator=(FilterDefinition&&) PBBAM_NOEXCEPT_MOVE_ASSIGN;
    ~FilterDefinition();

    const std::string& Id() const;
    const std::string& Description() const;

private:
    // required fields
    std::string id_;
    std::string description_;
};

///
/// \brief The FormatDefinition class
///
class FormatDefinition
{
public:
    FormatDefinition(std::string id, std::string number, std::string type, std::string description);

    FormatDefinition() = delete;
    FormatDefinition(const FormatDefinition&);
    FormatDefinition(FormatDefinition&&) noexcept;
    FormatDefinition& operator=(const FormatDefinition&);
    FormatDefinition& operator=(FormatDefinition&&) PBBAM_NOEXCEPT_MOVE_ASSIGN;
    ~FormatDefinition();

    const std::string& Id() const;
    const std::string& Number() const;
    const std::string& Type() const;
    const std::string& Description() const;

private:
    std::string id_;
    std::string number_;  // TODO: enum
    std::string type_;    // TODO: enum
    std::string description_;
};

///
/// \brief The GeneralDefinition class
///
class GeneralDefinition
{
public:
    GeneralDefinition(std::string id, std::string text);

    GeneralDefinition() = delete;
    GeneralDefinition(const GeneralDefinition&);
    GeneralDefinition(GeneralDefinition&&) noexcept;
    GeneralDefinition& operator=(const GeneralDefinition&);
    GeneralDefinition& operator=(GeneralDefinition&&) PBBAM_NOEXCEPT_MOVE_ASSIGN;
    ~GeneralDefinition();

    const std::string& Id() const;
    const std::string& Text() const;

private:
    // required fields
    std::string id_;
    std::string text_;
};

///
/// \brief The InfoDefinition class
///
class InfoDefinition
{
public:
    InfoDefinition(std::string id, std::string number, std::string type, std::string description,
                   std::string source = std::string{}, std::string version = std::string{});

    InfoDefinition() = delete;
    InfoDefinition(const InfoDefinition&);
    InfoDefinition(InfoDefinition&&) noexcept;
    InfoDefinition& operator=(const InfoDefinition&);
    InfoDefinition& operator=(InfoDefinition&&) PBBAM_NOEXCEPT_MOVE_ASSIGN;
    ~InfoDefinition();

    const std::string& Id() const;
    const std::string& Number() const;
    const std::string& Type() const;
    const std::string& Description() const;
    const boost::optional<std::string>& Source() const;
    const boost::optional<std::string>& Version() const;

    InfoDefinition& Source(std::string s);
    InfoDefinition& Version(std::string v);

private:
    // required fields
    // (functionally const, not marked as such to still allow moves)
    std::string id_;
    std::string number_;  // TODO: enum
    std::string type_;    // TODO: enum
    std::string description_;

    // optional fields - settable after ctor
    boost::optional<std::string> source_;
    boost::optional<std::string> version_;
};

}  // namespace VCF
}  // namespace PacBio

#endif  // PBBAM_VCF_VCFHEADERTYPES_H
