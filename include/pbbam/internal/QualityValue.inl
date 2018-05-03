// File Description
/// \file QualityValue.inl
/// \brief Inline implementations for the QualityValue class.
//
// Author: Derek Barnett

#include "pbbam/QualityValue.h"

namespace PacBio {
namespace BAM {

inline QualityValue::QualityValue(const uint8_t value)
    : value_{value}
{
    // clamp QV
    if (value_ > QualityValue::MAX)
        value_ = QualityValue::MAX;
}

inline char QualityValue::Fastq() const
{ return static_cast<char>(value_ + 33); }

inline QualityValue::operator uint8_t() const
{ return value_; }

inline QualityValue QualityValue::FromFastq(const char c)
{ return { static_cast<uint8_t>(c-33) }; }

} // namespace BAM
} // namespace PacBio
