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
/// \file QNameQuery.h
/// \brief Defines the QNameQuery class.
//
// Author: Derek Barnett

#ifndef QNAMEQUERY_H
#define QNAMEQUERY_H

#include "pbbam/internal/QueryBase.h"
#include <memory>

namespace PacBio {
namespace BAM {

/// \brief The QNameQuery class provides iterable access to a DataSet's records,
///        with each iteration of the query returning a contiguous block of
///        records that share a name.
///
/// There is no random-access here. It is simply a sequential read-through,
/// grouping contiguous results that share a BamRecord::FullName.
///
/// \note The name is not ideal - but for legacy reasons, it will remain as-is
///       for now. It will likely become something more explicit, like
///       "SequentialQNameGroupQuery", so that the name "QNameQuery" will be
///       available for a built-in query on a QNAME filter (or whitelist). This
///       will make it more consistent with other queries (ReadAccuracyQuery,
///       SubreadLengthQuery, ZmwQuery, etc).
///
class PBBAM_EXPORT QNameQuery : public internal::IGroupQuery
{
public:

    /// \brief Creates a new QNameQuery.
    ///
    /// \param[in] dataset      input data source(s)
    ///
    /// \throws std::runtime_error on failure to open/read underlying %BAM files
    ///
    QNameQuery(const DataSet& dataset);
    ~QNameQuery() override;

public:
    /// \brief Main iteration point for record access.
    ///
    /// Most client code should not need to use this method directly. Use
    /// iterators instead.
    ///
    bool GetNext(std::vector<BamRecord>& records) override;

private:
    struct QNameQueryPrivate;
    std::unique_ptr<QNameQueryPrivate> d_;
};

} // namespace BAM
} // namespace PacBio

#endif // QNAMEQUERY_H
