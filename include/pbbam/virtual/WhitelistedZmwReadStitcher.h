// Copyright (c) 2016, Pacific Biosciences of California, Inc.
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
/// \file WhitelistedZmwReadStitcher.h
/// \brief Defines the  ZmwReadStitcher class.
//
// Author: Derek Barnett

#ifndef WHITELISTEDZMWREADSTITCHER_H
#define WHITELISTEDZMWREADSTITCHER_H

#include "pbbam/Config.h"
#include "pbbam/virtual/VirtualZmwBamRecord.h"
#include <cstdint>
#include <memory>
#include <vector>
#include <string>

namespace PacBio {
namespace BAM {

class DataSet;
class PbiFilter;

/// \brief The WhitelistedZmwReadStitcher class provides an interface for
///        re-stitching "virtual" ZMW reads from their constituent parts,
///        limiting results to only those reads originating from a 'whitelist'
///         of ZMW hole numbers.
///
/// Whitelisted ZMWs that are not present in both primary and scraps BAMs
/// will be "pre-removed." This ensures that, given client code like this:
///
/// \include code/WhitelistedZmwReadStitcher.txt
///
/// each iteration will always provide valid data - either a valid virtual
/// record from Next() or a non-empty vector from NextRaw().
///
/// \note This reader requires that both input %BAM files also have associated
///       PBI files available for query. See BamFile::EnsurePacBioIndexExists .
///
class PBBAM_EXPORT WhitelistedZmwReadStitcher
{
public:
    /// \name Constructors & Related Methods
    /// \{

    /// \brief Creates a reader that will operate on a primary %BAM file (e.g. subread data)
    ///        and a scraps file, using a ZMW whitelist to filter the input.
    ///
    /// \param[in] zmwWhitelist         list of ZMWs to restrict iteration over
    /// \param[in] primaryBamFilePath   hqregion.bam or subreads.bam file path
    /// \param[in] scrapsBamFilePath    scraps.bam file path
    ///
    /// \note This reader requires that both input %BAM files also have associated PBI
    ///       files available for query. See BamFile::EnsurePacBioIndexExists .
    ///
    /// \throws std::runtime_error if any files (*.bam and/or *.pbi) were not available for reading, or
    ///         if malformed data encountered
    ///
    WhitelistedZmwReadStitcher(const std::vector<int32_t>& zmwWhitelist,
                              const std::string& primaryBamFilePath,
                              const std::string& scrapsBamFilePath);

    WhitelistedZmwReadStitcher() = delete;
    WhitelistedZmwReadStitcher(const WhitelistedZmwReadStitcher&) = delete;
    WhitelistedZmwReadStitcher(WhitelistedZmwReadStitcher&&)      = delete;
    WhitelistedZmwReadStitcher& operator=(const WhitelistedZmwReadStitcher&) = delete;
    WhitelistedZmwReadStitcher& operator=(WhitelistedZmwReadStitcher&&)      = delete;
    ~WhitelistedZmwReadStitcher();

    /// \}

public:
    /// \name Stitched Record Reading
    /// \{

    /// \returns true if more ZMWs are available for reading.
    bool HasNext() const;

    /// \returns the re-stitched polymerase read from the next ZMW in the whitelist
    VirtualZmwBamRecord Next();

    /// \returns the set of reads that belong to the next ZMW in the whitelist.
    ///          This enables stitching records in a distinct thread.
    ///
    std::vector<BamRecord> NextRaw();

    /// \}

public:
    /// \name File Headers
    /// \{

    /// \returns the BamHeader associated with this reader's "primary" %BAM file
    BamHeader PrimaryHeader() const;

    /// \returns the BamHeader associated with this reader's "scraps" %BAM file
    BamHeader ScrapsHeader() const;

    /// \}

private:
    struct WhitelistedZmwReadStitcherPrivate;
    std::unique_ptr<WhitelistedZmwReadStitcherPrivate> d_;
};

} // namespace BAM
} // namespace PacBio

#endif // WHITELISTEDZMWREADSTITCHER
