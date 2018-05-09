// File Description
/// \file Compare.cpp
/// \brief Implements the Compare class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/Compare.h"

#include <cstddef>
#include <functional>
#include <unordered_map>

namespace PacBio {
namespace BAM {
namespace internal {

struct TypeAlias
{
    std::string name_;
    std::string op_;
    std::string opAlpha_;

    TypeAlias(std::string name = std::string(), std::string op = std::string(),
              std::string opAlpha = std::string())
        : name_(std::move(name)), op_(std::move(op)), opAlpha_(std::move(opAlpha))
    {
    }
};

struct CompareTypeHash
{
    size_t operator()(const Compare::Type& t) const
    {
        return std::hash<int>()(static_cast<int>(t));
    }
};

// clang-format off
static const std::unordered_map<std::string, Compare::Type> opToTypeMap =
{
    // basic operators plus some permissiveness for other representations
    { "==",    Compare::EQUAL },
    { "=",     Compare::EQUAL },
    { "eq",    Compare::EQUAL },
    { "!=",    Compare::NOT_EQUAL },
    { "ne",    Compare::NOT_EQUAL },
    { "<",     Compare::LESS_THAN },
    { "lt",    Compare::LESS_THAN },
    { "&lt;",  Compare::LESS_THAN },
    { "<=",    Compare::LESS_THAN_EQUAL },
    { "lte",   Compare::LESS_THAN_EQUAL },
    { "&lt;=", Compare::LESS_THAN_EQUAL },
    { ">",     Compare::GREATER_THAN },
    { "gt",    Compare::GREATER_THAN },
    { "&gt;",  Compare::GREATER_THAN },
    { ">=",    Compare::GREATER_THAN_EQUAL },
    { "gte",   Compare::GREATER_THAN_EQUAL },
    { "&gt;=", Compare::GREATER_THAN_EQUAL },
    { "&",     Compare::CONTAINS },
    { "~",     Compare::NOT_CONTAINS }
};

static const std::unordered_map<Compare::Type, TypeAlias, CompareTypeHash> typeAliases =
{
    { Compare::EQUAL,              TypeAlias{ "Compare::EQUAL",              "==", "eq"  } },
    { Compare::NOT_EQUAL,          TypeAlias{ "Compare::NOT_EQUAL",          "!=", "ne"  } },
    { Compare::LESS_THAN,          TypeAlias{ "Compare::LESS_THAN",          "<",  "lt"  } },
    { Compare::LESS_THAN_EQUAL,    TypeAlias{ "Compare::LESS_THAN_EQUAL",    "<=", "lte" } },
    { Compare::GREATER_THAN,       TypeAlias{ "Compare::GREATER_THAN",       ">",  "gt"  } },
    { Compare::GREATER_THAN_EQUAL, TypeAlias{ "Compare::GREATER_THAN_EQUAL", ">=", "gte" } },
    { Compare::CONTAINS,           TypeAlias{ "Compare::CONTAINS",           "&",  "and" } },
    { Compare::NOT_CONTAINS,       TypeAlias{ "Compare::NOT_CONTAINS",       "~",  "not" } }
};
// clang-format on

}  // namespace internal

Compare::Type Compare::TypeFromOperator(const std::string& opString)
{
    try {
        return internal::opToTypeMap.at(opString);
    } catch (std::exception&) {
        throw std::runtime_error{opString + " is not a valid comparison operator."};
    }
}

std::string Compare::TypeToName(const Compare::Type& type)
{
    try {
        return internal::typeAliases.at(type).name_;
    } catch (std::exception&) {
        throw std::runtime_error{"invalid comparison type encountered"};
    }
}

std::string Compare::TypeToOperator(const Compare::Type& type, bool asAlpha)
{
    try {
        return asAlpha ? internal::typeAliases.at(type).opAlpha_
                       : internal::typeAliases.at(type).op_;
    } catch (std::exception&) {
        throw std::runtime_error{"invalid comparison type encountered"};
    }
}

}  // namespace BAM
}  // namespace PacBio
