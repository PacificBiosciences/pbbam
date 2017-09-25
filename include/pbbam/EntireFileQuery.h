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
/// \file EntireFileQuery.h
/// \brief Defines the EntireFileQuery class.
//
// Author: Derek Barnett

#ifndef ENTIREFILEQUERY_H
#define ENTIREFILEQUERY_H

#include "pbbam/internal/QueryBase.h"
#include <memory>

namespace PacBio {
namespace BAM {

/// \brief The EntireFileQuery class provides iterable access to a DataSet's
///        %BAM records, reading through the entire contents of each file.
///
/// Input files will be accessed in the order listed in the DataSet.
///
/// \include code/EntireFileQuery.txt
///
/// Iteration is not limited to only 'const' records. The files themselves will
/// not be affected, but individual records may be modified if needed.
///
/// \include code/EntireFileQuery_NonConst.txt
///
/// \note DataSets can be implicitly constructed from %BAM filenames as well.
///       Thus a single %BAM file can be read through using the following:
///
/// \include code/EntireFileQuery_BamFilename.txt
///
class PBBAM_EXPORT EntireFileQuery : public internal::IQuery
{
public:
    /// \brief Creates a new EntireFileQuery, reading through the entire
    ///        contents of a dataset.
    ///
    /// \param[in] dataset  input data source(s)
    /// \throws std::runtime_error on failure to open/read underlying %BAM
    ///         files.
    ///
    EntireFileQuery(const PacBio::BAM::DataSet& dataset);
    ~EntireFileQuery() override;

public:
    /// \brief Main iteration point for record access.
    ///
    /// Most client code should not need to use this method directly. Use
    /// iterators instead.
    ///
    bool GetNext(BamRecord& r) override;

private:
    struct EntireFileQueryPrivate;
    std::unique_ptr<EntireFileQueryPrivate> d_;
};

} // namespace BAM
} // namspace PacBio

#endif // ENTIREFILEQUERY_H
