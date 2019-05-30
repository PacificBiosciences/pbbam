// File Description
/// \file QualityValue.h
/// \brief Implements the QualityValue class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/QualityValue.h"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <type_traits>

namespace PacBio {
namespace BAM {

const uint8_t QualityValue::MAX = 93;

static_assert(std::is_copy_constructible<QualityValue>::value,
              "QualityValue(const QualityValue&) is not = default");
static_assert(std::is_copy_assignable<QualityValue>::value,
              "QualityValue& operator=(const QualityValue&) is not = default");

static_assert(std::is_nothrow_move_constructible<QualityValue>::value,
              "QualityValue(QualityValue&&) is not = noexcept");
static_assert(std::is_nothrow_move_assignable<QualityValue>::value,
              "QualityValue& operator=(QualityValue&&) is not = noexcept");

QualityValue::QualityValue(const uint8_t value)
    : value_{// clamp QV
             std::min(value, QualityValue::MAX)}
{
}

char QualityValue::Fastq() const { return static_cast<char>(value_ + 33); }

QualityValue::operator uint8_t() const { return value_; }

QualityValue QualityValue::FromFastq(const char c) { return {static_cast<uint8_t>(c - 33)}; }

}  // namespace BAM
}  // namespace PacBio
