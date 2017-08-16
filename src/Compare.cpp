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
//
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

    TypeAlias(const std::string& name = std::string(),
              const std::string& op = std::string(),
              const std::string& opAlpha = std::string())
        : name_(name)
        , op_(op)
        , opAlpha_(opAlpha)
    { }
};

struct CompareTypeHash
{
    size_t operator()(const Compare::Type& t) const
    { return std::hash<int>()(static_cast<int>(t)); }
};

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

} // namespace internal

Compare::Type Compare::TypeFromOperator(const std::string& opString)
{
    try {
        return internal::opToTypeMap.at(opString);
    } catch (std::exception&) {
        throw std::runtime_error(opString + " is not a valid comparison operator." );
    }
}

std::string Compare::TypeToName(const Compare::Type& type)
{
    try {
        return internal::typeAliases.at(type).name_;
    } catch (std::exception&) {
        throw std::runtime_error("invalid comparison type encountered" );
    }
}

std::string Compare::TypeToOperator(const Compare::Type& type, bool asAlpha)
{
    try {
        return asAlpha ? internal::typeAliases.at(type).opAlpha_
                       : internal::typeAliases.at(type).op_;
    } catch (std::exception&) {
        throw std::runtime_error("invalid comparison type encountered" );
    }
}

} // namespace BAM
} // namespace PacBio
