// File Description
/// \file QNameQuery.cpp
/// \brief Implements the QNameQuery class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/QualityValues.h"

#include <algorithm>
#include <cassert>
#include <type_traits>

#include <boost/algorithm/string.hpp>

namespace PacBio {
namespace BAM {

static_assert(std::is_copy_constructible<QualityValues>::value,
              "QualityValues(const QualityValues&) is not = default");
static_assert(std::is_copy_assignable<QualityValues>::value,
              "QualityValues& operator=(const QualityValues&) is not = default");

static_assert(std::is_nothrow_move_constructible<QualityValues>::value,
              "QualityValues(QualityValues&&) is not = noexcept");
static_assert(std::is_nothrow_move_assignable<QualityValues>::value,
              "QualityValues& operator=(QualityValues&&) is not = noexcept");

QualityValues::QualityValues(const std::string& fastqString) : std::vector<QualityValue>{}
{
    std::string fastqString_{std::move(fastqString)};
    boost::algorithm::trim(fastqString_);
    resize(fastqString_.size());
    std::transform(fastqString_.cbegin(), fastqString_.cend(), begin(), QualityValue::FromFastq);
}

QualityValues::QualityValues(std::vector<QualityValue> quals)
    : std::vector<QualityValue>{std::move(quals)}
{
}

QualityValues::QualityValues(const std::vector<uint8_t>& quals)
    : std::vector<QualityValue>(quals.size())
{
    std::copy(quals.cbegin(), quals.cend(), begin());
}

QualityValues::QualityValues(const std::vector<uint8_t>::const_iterator first,
                             const std::vector<uint8_t>::const_iterator last)
    : std::vector<QualityValue>(first, last)
{
}

QualityValues::QualityValues(const QualityValues::const_iterator first,
                             const QualityValues::const_iterator last)
    : std::vector<QualityValue>{}
{
    assign(first, last);
}

QualityValues& QualityValues::operator=(std::vector<QualityValue> quals)
{
    std::vector<QualityValue>::operator=(std::move(quals));
    return *this;
}

std::vector<QualityValue>::const_iterator QualityValues::cbegin() const
{
    return std::vector<QualityValue>::cbegin();
}

std::vector<QualityValue>::const_iterator QualityValues::cend() const
{
    return std::vector<QualityValue>::cend();
}

std::vector<QualityValue>::const_iterator QualityValues::begin() const
{
    return std::vector<QualityValue>::begin();
}

std::vector<QualityValue>::const_iterator QualityValues::end() const
{
    return std::vector<QualityValue>::end();
}

std::vector<QualityValue>::iterator QualityValues::begin()
{
    return std::vector<QualityValue>::begin();
}

std::vector<QualityValue>::iterator QualityValues::end()
{
    return std::vector<QualityValue>::end();
}

QualityValues QualityValues::FromFastq(const std::string& fastq) { return QualityValues{fastq}; }

std::string QualityValues::Fastq() const
{
    std::string result;
    result.reserve(size());
    for (const auto qv : *this)
        result.push_back(qv.Fastq());
    return result;
}

bool QualityValues::operator==(const std::string& fastq) const
{
    return *this == QualityValues(fastq);
}

bool QualityValues::operator!=(const std::string& fastq) const
{
    return *this != QualityValues(fastq);
}

}  // namespace BAM
}  // namespace PacBio
