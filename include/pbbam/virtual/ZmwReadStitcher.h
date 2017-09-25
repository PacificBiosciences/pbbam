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
/// \file ZmwReadStitcher.h
/// \brief Defines the ZmwReadStitcher class.
//
// Author: Derek Barnett

#ifndef ZMWREADSTITCHER_H
#define ZMWREADSTITCHER_H

#include "pbbam/BamRecord.h"
#include "pbbam/Config.h"
#include "pbbam/virtual/VirtualZmwBamRecord.h"
#include <memory>
#include <vector>
#include <string>

namespace PacBio {
namespace BAM {

class DataSet;
class PbiFilter;

/// \brief The ZmwReadStitcher class provides an interface for re-stitching
///        "virtual" polymerase reads from their constituent parts.
///
/// \note This reader requires that any input %BAM files also have associated PBI
///       files available for query. See BamFile::EnsurePacBioIndexExists .
///
class PBBAM_EXPORT ZmwReadStitcher
{
public:
    /// \name Constructors & Related Methods
    /// \{

    /// entire file, from BAM names
    ZmwReadStitcher(const std::string& primaryBamFilePath,
                    const std::string& scrapsBamFilePath);

    /// filtered input from BAM names
    ZmwReadStitcher(const std::string& primaryBamFilePath,
                    const std::string& scrapsBamFilePath,
                    const PbiFilter& filter);

    /// maybe filtered, from DataSet input
    ZmwReadStitcher(const DataSet& dataset);

    ZmwReadStitcher() = delete;
    ZmwReadStitcher(const ZmwReadStitcher&) = delete;
    ZmwReadStitcher(ZmwReadStitcher&&) = delete;
    ZmwReadStitcher& operator=(const ZmwReadStitcher&) = delete;
    ZmwReadStitcher& operator=(ZmwReadStitcher&&) = delete;
    ~ZmwReadStitcher();

    /// \}

public:
    /// \name File Headers
    /// \{

    /// \returns the BamHeader associated with this reader's "primary" %BAM file
    BamHeader PrimaryHeader() const;

    /// \returns the BamHeader associated with this reader's "scraps" %BAM file
    BamHeader ScrapsHeader() const;

    /// \}

public:
    /// \name Stitched Record Reading
    ///

    /// \returns true if more ZMWs are available for reading.
    bool HasNext();

    /// \returns the next stitched polymerase read
    VirtualZmwBamRecord Next();

    /// \returns the next set of reads that belong to one ZMW.
    ///          This enables stitching records in a distinct thread.
    ///
    std::vector<BamRecord> NextRaw();

    /// \}

private:
    struct ZmwReadStitcherPrivate;
    std::unique_ptr<ZmwReadStitcherPrivate> d_;
};

} // namespace BAM
} // namespace PacBio

#endif // ZMWREADSTITCHER_H
