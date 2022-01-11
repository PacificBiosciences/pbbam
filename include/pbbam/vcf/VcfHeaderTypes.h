#ifndef PBBAM_VCF_VCFHEADERTYPES_H
#define PBBAM_VCF_VCFHEADERTYPES_H

#include <pbbam/Config.h>

#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace PacBio {
namespace VCF {

using Sample = std::string;

class ContigDefinition
{
public:
    explicit ContigDefinition(std::string id);
    ContigDefinition(std::string id, std::vector<std::pair<std::string, std::string>> attributes);

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

    const std::string& Id() const;
    const std::string& Number() const;
    const std::string& Type() const;
    const std::string& Description() const;
    const std::optional<std::string>& Source() const;
    const std::optional<std::string>& Version() const;

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
    std::optional<std::string> source_;
    std::optional<std::string> version_;
};

}  // namespace VCF
}  // namespace PacBio

#endif  // PBBAM_VCF_VCFHEADERTYPES_H
