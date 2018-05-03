// File Description
/// \file QualityValues.inl
/// \brief Inline implementations for the QualityValues class.
//
// Author: Derek Barnett

#include "pbbam/QualityValues.h"
#include <algorithm>

namespace PacBio {
namespace BAM {

inline QualityValues::QualityValues(const std::string& fastqString)
    : std::vector<QualityValue>{}
{
    resize(fastqString.size());
    std::transform(fastqString.cbegin(), fastqString.cend(),
                   begin(), QualityValue::FromFastq);
}

inline QualityValues::QualityValues(std::vector<QualityValue> quals)
    : std::vector<QualityValue>{std::move(quals)}
{ }

inline QualityValues::QualityValues(const std::vector<uint8_t>& quals)
    : std::vector<QualityValue>(quals.size())
{
    std::copy(quals.cbegin(), quals.cend(), begin());
}

inline QualityValues::QualityValues(const std::vector<uint8_t>::const_iterator first,
                                    const std::vector<uint8_t>::const_iterator last)
    : std::vector<QualityValue>(first, last)
{ }

inline QualityValues::QualityValues(const QualityValues::const_iterator first,
                                    const QualityValues::const_iterator last)
    : std::vector<QualityValue>{}
{
    assign(first, last);
}

inline QualityValues& QualityValues::operator=(std::vector<QualityValue> quals)
{ std::vector<QualityValue>::operator=(std::move(quals)); return *this; }

inline std::vector<QualityValue>::const_iterator QualityValues::cbegin() const
{ return std::vector<QualityValue>::cbegin(); }

inline std::vector<QualityValue>::const_iterator QualityValues::cend() const
{ return std::vector<QualityValue>::cend(); }

inline std::vector<QualityValue>::const_iterator QualityValues::begin() const
{ return std::vector<QualityValue>::begin(); }

inline std::vector<QualityValue>::const_iterator QualityValues::end() const
{ return std::vector<QualityValue>::end(); }

inline std::vector<QualityValue>::iterator QualityValues::begin()
{ return std::vector<QualityValue>::begin(); }

inline std::vector<QualityValue>::iterator QualityValues::end()
{ return std::vector<QualityValue>::end(); }

inline QualityValues QualityValues::FromFastq(const std::string& fastq)
{ return QualityValues{fastq}; }

inline std::string QualityValues::Fastq() const
{
    std::string result;
    result.reserve(size());
    for (const auto qv : *this)
        result.push_back(qv.Fastq());
    return result;
}

inline bool QualityValues::operator==(const std::string& fastq) const
{ return *this == QualityValues(fastq); }

inline bool QualityValues::operator!=(const std::string& fastq) const
{ return *this != QualityValues(fastq); }

} // namespace BAM
} // namespace PacBio
