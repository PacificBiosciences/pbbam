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

// Author: Derek Barnett

#ifndef ZMWWHITELISTVIRTUALREADER_H
#define ZMWWHITELISTVIRTUALREADER_H

#include <deque>
#include <memory>
#include <vector>
#include <string>
#include "pbbam/BamFile.h"
#include "pbbam/BamRecord.h"
#include "pbbam/Config.h"
#include "pbbam/PbiIndexedBamReader.h"
#include "pbbam/virtual/VirtualPolymeraseBamRecord.h"

namespace PacBio {
namespace BAM {

class ZmwWhitelistVirtualReader
{
public:
    /// \name Constructors & Related Methods
    /// \{

    /// Constructor takes a zmw whitelist & two input BAM file paths.
    ///
    /// Whitelisted ZMWs that are not present in both primary and scraps BAMs
    /// will be "pre-removed." This ensures that given client code like this:
    ///
    /// \code{.cpp}
    /// ZmwWhitelistVirtualReader reader(zmws, f1, f2);
    /// while(reader.HasNext()) {
    ///     auto virtualRecord = reader.Next();
    ///     // do stuf...
    /// }
    /// \endcode
    ///
    /// each iteration will always provide valid data: either a valid virtual record from
    /// Next() or a non-empty vector from NextRaw().
    ///
    /// \note This reader requires that both input BAM files also have associated PBI
    ///       files available for query. See BamFile::EnsurePacBioIndexExists .
    ///
    /// \param[in] zmwWhitelist       list of ZMWs to restrict iteration over
    /// \param[in] primaryBamFilePath hqregion.bam or subreads.bam file path
    /// \param[in] scrapsBamFilePath  scraps.bam file path
    ///
    /// \throws std::runtime_error if any files (*.bam and/or *.pbi) were not available for reading, or
    ///         if malformed data encountered
    ///
    ZmwWhitelistVirtualReader(const std::vector<int32_t>& zmwWhitelist,
                              const std::string& primaryBamFilePath,
                              const std::string& scrapsBamFilePath);

    ZmwWhitelistVirtualReader(void) = delete;
    ZmwWhitelistVirtualReader(const ZmwWhitelistVirtualReader&) = delete;
    ZmwWhitelistVirtualReader(ZmwWhitelistVirtualReader&&)      = delete;
    ZmwWhitelistVirtualReader& operator=(const ZmwWhitelistVirtualReader&) = delete;
    ZmwWhitelistVirtualReader& operator=(ZmwWhitelistVirtualReader&&)      = delete;
    ~ZmwWhitelistVirtualReader(void) = default;

    /// \}

public:
    /// \name Stitched Record Reading
    /// \{

    /// Provides the stitched polymerase read from the next ZMW (with data) in the whitelist
    VirtualPolymeraseBamRecord Next(void);

    /// Provides the set of reads that belong to next ZMW (with data) in the whitelist.
    /// Enables stitching records in a distinct thread.
    std::vector<BamRecord> NextRaw(void);

    /// \returns true if more ZMWs are available for reading.
    bool HasNext(void) const;

    /// \}

public:
    /// \name File Metadata
    /// \{

    BamHeader PrimaryHeader(void) const;
    BamHeader ScrapsHeader(void) const;

    /// \}

private:
    const std::string        primaryBamFilePath_;
    const std::string        scrapsBamFilePath_;
    std::unique_ptr<BamFile> primaryBamFile_;
    std::unique_ptr<BamFile> scrapsBamFile_;
    std::unique_ptr<PbiIndexedBamReader> primaryReader_;
    std::unique_ptr<PbiIndexedBamReader> scrapsReader_;
    std::unique_ptr<BamHeader> polyHeader_;
    std::deque<int32_t>        zmwWhitelist_;
};

} // namespace BAM
} // namespace PacBio

#endif // ZMWWHITELISTVIRTUALREADER_H
