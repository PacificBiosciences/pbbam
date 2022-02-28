#include "PbbamInternalConfig.h"

#include <pbbam/Compare.h>

#include <functional>
#include <unordered_map>

#include <cstddef>

namespace PacBio {
namespace BAM {
namespace {

struct TypeAlias
{
    std::string name_;
    std::string op_;
    std::string opAlpha_;

    TypeAlias(std::string name = std::string(), std::string op = std::string(),
              std::string opAlpha = std::string())
        : name_(std::move(name)), op_(std::move(op)), opAlpha_(std::move(opAlpha))
    {}
};

struct CompareTypeHash
{
    size_t operator()(const Compare::Type& t) const noexcept
    {
        return std::hash<int>()(static_cast<int>(t));
    }
};

// clang-format off
static const std::unordered_map<std::string, Compare::Type> opToTypeMap =
{
    // basic operators plus some permissiveness for other representations
    { "==",     Compare::EQUAL },
    { "=",      Compare::EQUAL },
    { "eq",     Compare::EQUAL },
    { "in",     Compare::EQUAL },
    { "!=",     Compare::NOT_EQUAL },
    { "ne",     Compare::NOT_EQUAL },
    { "not_in", Compare::NOT_EQUAL },
    { "<",      Compare::LESS_THAN },
    { "lt",     Compare::LESS_THAN },
    { "&lt;",   Compare::LESS_THAN },
    { "<=",     Compare::LESS_THAN_EQUAL },
    { "lte",    Compare::LESS_THAN_EQUAL },
    { "&lt;=",  Compare::LESS_THAN_EQUAL },
    { ">",      Compare::GREATER_THAN },
    { "gt",     Compare::GREATER_THAN },
    { "&gt;",   Compare::GREATER_THAN },
    { ">=",     Compare::GREATER_THAN_EQUAL },
    { "gte",    Compare::GREATER_THAN_EQUAL },
    { "&gt;=",  Compare::GREATER_THAN_EQUAL },
    { "&",      Compare::CONTAINS },
    { "~",      Compare::NOT_CONTAINS }
};

static const std::unordered_map<Compare::Type, TypeAlias, CompareTypeHash> typeAliases =
{
    { Compare::EQUAL,              TypeAlias{ "Compare::EQUAL",              "==", "eq" } },
    { Compare::NOT_EQUAL,          TypeAlias{ "Compare::NOT_EQUAL",          "!=", "ne" } },
    { Compare::LESS_THAN,          TypeAlias{ "Compare::LESS_THAN",          "<",  "lt"  } },
    { Compare::LESS_THAN_EQUAL,    TypeAlias{ "Compare::LESS_THAN_EQUAL",    "<=", "lte" } },
    { Compare::GREATER_THAN,       TypeAlias{ "Compare::GREATER_THAN",       ">",  "gt"  } },
    { Compare::GREATER_THAN_EQUAL, TypeAlias{ "Compare::GREATER_THAN_EQUAL", ">=", "gte" } },
    { Compare::CONTAINS,           TypeAlias{ "Compare::CONTAINS",           "&",  "and" } },
    { Compare::NOT_CONTAINS,       TypeAlias{ "Compare::NOT_CONTAINS",       "~",  "not" } }
};
// clang-format on

}  // namespace

bool Compare::AlignmentPosition::operator()(const BamRecord& lhs, const BamRecord& rhs) const
{
    const int32_t lhsId = lhs.ReferenceId();
    const int32_t rhsId = rhs.ReferenceId();

    // push unmapped reads to bottom
    if (lhsId == -1) {
        return false;
    }
    if (rhsId == -1) {
        return true;
    }

    // compare by refId, then position
    if (lhsId == rhsId) {
        return lhs.ReferenceStart() < rhs.ReferenceStart();
    } else {
        return lhsId < rhsId;
    }
}

bool Compare::QName::operator()(const BamRecord& lhs, const BamRecord& rhs) const
{
    // movie name
    const auto lMovieName = lhs.MovieName();
    const auto rMovieName = rhs.MovieName();
    const int cmp = lMovieName.compare(rMovieName);
    if (cmp != 0) {
        return cmp < 0;
    }

    // hole number
    const auto lhsZmw = lhs.HoleNumber();
    const auto rhsZmw = rhs.HoleNumber();
    if (lhsZmw != rhsZmw) {
        return lhsZmw < rhsZmw;
    }

    // shuffle CCS/transcript reads after all others
    if (IsCcsOrTranscript(lhs.Type())) {
        return false;
    }
    if (IsCcsOrTranscript(rhs.Type())) {
        return true;
    }

    // sort on qStart, then finally qEnd
    const auto lhsQStart = lhs.QueryStart();
    const auto rhsQStart = rhs.QueryStart();
    if (lhsQStart != rhsQStart) {
        return lhsQStart < rhsQStart;
    }

    const auto lhsQEnd = lhs.QueryEnd();
    const auto rhsQEnd = rhs.QueryEnd();
    return lhsQEnd < rhsQEnd;
}

Compare::Type Compare::TypeFromOperator(const std::string& opString)
{
    try {
        return opToTypeMap.at(opString);
    } catch (std::exception&) {
        throw std::runtime_error{"[pbbam] comparison ERROR: " + opString +
                                 " is not a valid comparison operator."};
    }
}

std::string Compare::TypeToName(const Compare::Type& type)
{
    try {
        return typeAliases.at(type).name_;
    } catch (std::exception&) {
        throw std::runtime_error{"[pbbam] comparison ERROR: invalid comparison type encountered"};
    }
}

std::string Compare::TypeToOperator(const Compare::Type& type, bool asAlpha)
{
    try {
        return asAlpha ? typeAliases.at(type).opAlpha_ : typeAliases.at(type).op_;
    } catch (std::exception&) {
        throw std::runtime_error{"[pbbam] comparison ERROR: invalid comparison type encountered"};
    }
}

}  // namespace BAM
}  // namespace PacBio
