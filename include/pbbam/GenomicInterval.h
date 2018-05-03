// File Description
/// \file GenomicInterval.h
/// \brief Defines the GenomicInterval class.
//
// Author: Derek Barnett

#ifndef GENOMICINTERVAL_H
#define GENOMICINTERVAL_H

#include <cstddef>
#include <string>
#include "pbbam/Config.h"
#include "pbbam/Interval.h"
#include "pbbam/Position.h"

namespace PacBio {
namespace BAM {

/// \brief The GenomicInterval class represents a genomic interval (reference
///        name and 0-based coordinates).
///
class PBBAM_EXPORT GenomicInterval
{
public:
    /// \name Constructors & Related Methods
    ///  \{

    /// \brief Creates an empty genomic interval
    GenomicInterval() = default;

    /// \brief Creates a genomic interval on sequence with \p name, using range:
    ///       [\p start, \p stop)
    GenomicInterval(std::string name, Position start, Position stop);

    /// \brief Creates a genomic interval, using REGION string
    ///
    /// "<ref>:<start>-<stop>" ("chr8:200-600")
    ///
    /// \note The htslib/samtools REGION string expects start positions to be
    ///       1-based. However, throughout pbbam (including the rest of this
    ///       class), we stick to 0-based start coordinates. Thus, while the
    ///       syntax matches that of samtools, we are using a 0-based start
    ///       coordinate here.
    ///
    GenomicInterval(const std::string& zeroBasedRegionString);

    GenomicInterval(const GenomicInterval&) = default;
    GenomicInterval(GenomicInterval&&) = default;
    GenomicInterval& operator=(const GenomicInterval&) = default;
    GenomicInterval& operator=(GenomicInterval&&) = default;
    ~GenomicInterval() = default;

    /// \}

public:
    /// \name Comparison Operators
    /// \{

    /// \returns true if same id & underlying interval
    bool operator==(const GenomicInterval& other) const;

    /// \returns true if either ids or underlying intervals differ
    bool operator!=(const GenomicInterval& other) const;

    /// \}

public:
    /// \name Interval Operations
    /// \{

    /// \returns true if same id and underlying Interval::CoveredBy() other.
    bool CoveredBy(const GenomicInterval& other) const;

    /// \returns true if same id and underlying Interval::Covers() other.
    bool Covers(const GenomicInterval& other) const;

    /// \returns true if same id and underlying Interval::Intersects() other.
    bool Intersects(const GenomicInterval& other) const;

    /// \returns true if underlying Interval::IsValid(), and id/endpoints are
    ///          non-negative.
    ///
    bool IsValid() const;

    /// \returns length of underlying
    size_t Length() const;

    /// \}

public:
    /// \name Attributes
    /// \{

    /// \returns interval reference name
    std::string Name() const;

    /// \returns underlying Interval object
    PacBio::BAM::Interval<Position> Interval() const;

    /// \returns interval start coordinate
    Position Start() const;

    /// \returns interval stop coordinate
    Position Stop() const;

    /// \}

public:
    /// \name Attributes
    /// \{

    /// Sets this interval's reference name.
    ///
    /// \param[in] name
    /// \returns reference to this interval
    ///
    GenomicInterval& Name(std::string name);

    /// Sets this underlying Interval
    ///
    /// \param[in] interval
    /// \returns reference to this interval
    ///
    GenomicInterval& Interval(PacBio::BAM::Interval<Position> interval);

    /// Sets this interval's start coordinate.
    ///
    /// \param[in] start
    /// \returns reference to this interval
    ///
    GenomicInterval& Start(const Position start);

    /// Sets this interval's stop coordinate.
    ///
    /// \param[in] stop
    /// \returns reference to this interval
    ///
    GenomicInterval& Stop(const Position stop);

    /// \}

private:
    std::string name_;
    PacBio::BAM::Interval<Position> interval_;
};

}  // namespace BAM
}  // namspace PacBio

#include "pbbam/internal/GenomicInterval.inl"

#endif  // GENOMICINTERVAL_H
