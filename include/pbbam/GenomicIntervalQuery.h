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
/// \file GenomicIntervalQuery.h
/// \brief Defines the GenomicIntervalQuery class.
//
// Author: Derek Barnett

#ifndef GENOMICINTERVALQUERY_H
#define GENOMICINTERVALQUERY_H

#include "pbbam/GenomicInterval.h"
#include "pbbam/internal/QueryBase.h"
#include <memory>

namespace PacBio {
namespace BAM {

/// \brief The GenomicIntervalQuery class provides iterable access to a
///        DataSet's %BAM records, limiting results to those overlapping a
///        GenomicInterval.
///
/// Example:
/// \include code/GenomicIntervalQuery.txt
///
/// \note Currently, all %BAM files must have a corresponding ".bai" index file.
///       Use BamFile::EnsureStandardIndexExists before creating the query if
///       one may not be present.
///
class PBBAM_EXPORT GenomicIntervalQuery : public internal::IQuery
{
public:

    /// \brief Constructs a new GenomiIntervalQuery, limiting record results to
    ///        only those overalpping a GenomicInterval.
    ///
    /// \param[in] interval genomic interval of interest
    /// \param[in] dataset  input data source(s)
    ///
    /// \throws std::runtime_error on failure to open/read underlying %BAM or
    ///         BAI files.
    ///
    GenomicIntervalQuery(const GenomicInterval& interval,
                         const PacBio::BAM::DataSet& dataset);
    ~GenomicIntervalQuery() override;

public:
    /// \brief Main iteration point for record access.
    ///
    /// Most client code should not need to use this method directly. Use
    /// iterators instead.
    ///
    bool GetNext(BamRecord& r) override;

public:
    /// \brief Sets a new genomic interval.
    ///
    /// This allows the same dataset/query to be re-used over multiple regions of
    /// interest:
    ///
    /// \include code/GenomicIntervalQuery_Reuse.txt
    ///
    /// \param[in] interval new genomic interval
    /// \returns reference to this query
    ///
    GenomicIntervalQuery& Interval(const GenomicInterval& interval);

    /// \returns Current genomic interval active on this query.
    const GenomicInterval& Interval() const;

private:
    struct GenomicIntervalQueryPrivate;
    std::unique_ptr<GenomicIntervalQueryPrivate> d_;
};

} // namespace BAM
} // namspace PacBio

#endif // GENOMICINTERVALQUERY_H
