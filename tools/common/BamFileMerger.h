// Copyright (c) 2015, Pacific Biosciences of California, Inc.
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
// DISCLAIMED. IN NO EVENT  SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

// Author: Derek Barnett

#ifndef BAMFILEMERGER_H
#define BAMFILEMERGER_H

#include <pbbam/DataSet.h>
#include <pbbam/PbiFilter.h>
#include <pbbam/ProgramInfo.h>
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {
namespace common {

class BamFileMerger
{
public:
    /// \brief Runs merger on a dataset, applying any supplied filters.
    ///
    /// When this function exits, a merged BAM (and optional PBI) will have been
    /// written and closed.
    ///
    /// \param[in] dataset          provides input filenames & filters
    /// \param[in] outputFilename   resulting BAM output
    /// \param[in] mergeProgram     info about the calling program. Adds a @PG entry to merged header.
    /// \param[in] createPbi        if true, creates a PBI alongside output BAM
    ///
    /// \throws std::runtime_error if any any errors encountered while reading or writing
    ///
    static void Merge(const PacBio::BAM::DataSet& dataset,
                      const std::string& outputFilename,
                      const PacBio::BAM::ProgramInfo& mergeProgram = PacBio::BAM::ProgramInfo(),
                      bool createPbi = true);
};

} // namespace common
} // namespace BAM
} // namespace PacBio

#include "BamFileMerger.inl"

#endif // BAMFILEMERGER_H
