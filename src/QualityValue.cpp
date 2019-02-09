// File Description
/// \file QualityValue.h
/// \brief Implements the QualityValue class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/QualityValue.h"

#include <cstdint>

namespace PacBio {
namespace BAM {

const uint8_t QualityValue::MAX = 93;

QualityValue::QualityValue(const uint8_t value) : value_{value}
{
    // clamp QV
    if (value_ > QualityValue::MAX) value_ = QualityValue::MAX;
}

char QualityValue::Fastq() const { return static_cast<char>(value_ + 33); }

QualityValue::operator uint8_t() const { return value_; }

QualityValue QualityValue::FromFastq(const char c) { return {static_cast<uint8_t>(c - 33)}; }

}  // namespace BAM
}  // namespace PacBio
