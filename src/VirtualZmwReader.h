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
/// \file VirtualZmwReader.h
/// \brief Defines the VirtualZmwReader class.
//
// Author: Armin Töpfer

#ifndef VIRTUALZMWREADER_H
#define VIRTUALZMWREADER_H

#include <memory>

#include "pbbam/BamFile.h"
#include "pbbam/BamRecord.h"
#include "pbbam/Config.h"
#include "pbbam/EntireFileQuery.h"
#include "pbbam/PbiFilter.h"
#include "pbbam/PbiFilterQuery.h"
#include "pbbam/virtual/VirtualZmwBamRecord.h"

namespace PacBio {
namespace BAM {
namespace internal {

class VirtualZmwReader
{
public:
    /// \brief Creates a reader that will operate on a primary %BAM file (e.g.
    ///        subread data) and a scraps file, consuming all reads.
    ///
    /// \param[in] primaryBamFilepath hqregion.bam or subreads.bam file path
    /// \param[in] scrapsBamFilepath  scraps.bam file path
    ///
    VirtualZmwReader(const std::string& primaryBamFilepath,
                     const std::string& scrapsBamFilepath);

    /// \brief Creates a reader that will operate on a primary %BAM file (e.g.
    ///        subread data) and a scraps file, respecting the provided PBI
    ///        filter.
    ///
    /// \note All %BAM files must have a corresponding ".pbi" index file to use
    ///       the filter. You may need to call BamFile::EnsurePacBioIndexExists
    ///       before constructing the reader.
    ///
    /// \param[in] primaryBamFilepath hqregion.bam or subreads.bam file path
    /// \param[in] scrapsBamFilepath  scraps.bam file path
    /// \param[in] filter PBI filter criteria
    ///
    VirtualZmwReader(const std::string& primaryBamFilepath,
                     const std::string& scrapsBamFilepath,
                     const PbiFilter& filter);

    VirtualZmwReader(void) = delete;
    VirtualZmwReader(const VirtualZmwReader&) = delete;
    VirtualZmwReader(VirtualZmwReader&&) = delete;
    VirtualZmwReader& operator=(const VirtualZmwReader&) = delete;
    VirtualZmwReader& operator=(VirtualZmwReader&&) = delete;
    ~VirtualZmwReader(void);

public:

    /// \returns the BamHeader associated with this reader's "primary" %BAM file
    BamHeader PrimaryHeader(void) const;

    /// \returns the BamHeader associated with this reader's "scraps" %BAM file
    BamHeader ScrapsHeader(void) const;

public:

    /// \returns true if more ZMWs are available for reading.
    bool HasNext(void);

    /// \returns the next stitched polymerase read
    VirtualZmwBamRecord Next(void);

    /// \returns the next set of reads that belong to one ZMW.
    ///          This enables stitching records in a distinct thread.
    ///
    std::vector<BamRecord> NextRaw(void);

private:
    std::unique_ptr<BamFile>          primaryBamFile_;
    std::unique_ptr<BamFile>          scrapsBamFile_;
    std::unique_ptr<internal::IQuery> primaryQuery_;
    std::unique_ptr<internal::IQuery> scrapsQuery_;
    internal::IQuery::iterator        primaryIt_;
    internal::IQuery::iterator        scrapsIt_;
    std::unique_ptr<BamHeader>        stitchedHeader_;
};

} // namespace internal
} // namespace BAM
} // namespace PacBio

#endif // VirtualZmwREADER_H
