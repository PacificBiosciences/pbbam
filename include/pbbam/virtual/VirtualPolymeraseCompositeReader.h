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
/// \file VirtualPolymeraseCompositeReader.h
/// \brief Defines the VirtualPolymeraseCompositeReader class.
//
// Author: Derek Barnett

#ifndef VIRTUALPOLYMERASECOMPOSITEREADER_H
#define VIRTUALPOLYMERASECOMPOSITEREADER_H

#include "pbbam/DataSet.h"
#include "pbbam/virtual/VirtualPolymeraseReader.h"
#include <deque>
#include <memory>
#include <string>
#include <utility>

namespace PacBio {
namespace BAM {

/// \brief The VirtualPolymeraseCompositeReader provides an interface for
///        re-stitching "virtual" polymerase reads from their constituent parts,
///        across multiple %BAM resources from a DataSet.
///
/// This class is essentially a DataSet-aware wrapper around
/// VirtualPolymeraseReader, enabling multiple resources as input. See that
/// class's documentation for more info.
///
class PBBAM_EXPORT VirtualPolymeraseCompositeReader
{
public:
    /// \name Constructors & Related Methods
    /// \{

    VirtualPolymeraseCompositeReader(const DataSet& dataset);

    VirtualPolymeraseCompositeReader(void) = delete;
    VirtualPolymeraseCompositeReader(const VirtualPolymeraseCompositeReader&) = delete;
    VirtualPolymeraseCompositeReader(VirtualPolymeraseCompositeReader&&) = delete;
    VirtualPolymeraseCompositeReader& operator=(const VirtualPolymeraseCompositeReader&) = delete;
    VirtualPolymeraseCompositeReader& operator=(VirtualPolymeraseCompositeReader&&) = delete;
    ~VirtualPolymeraseCompositeReader(void) = default;

    /// \}

public:
    /// \name Stitched Record Reading
    ///

    /// \returns true if more ZMWs/files are available for reading.
    bool HasNext(void);

    /// \returns the next stitched polymerase read
    VirtualPolymeraseBamRecord Next(void);

    /// \returns the next set of reads that belong to one ZMW from one %BAM
    ///          resource (a primary %BAM and/or its scraps file). This enables
    ///          stitching records in a distinct thread.
    ///
    std::vector<BamRecord> NextRaw(void);

    /// \}

private:
    std::deque< std::pair<std::string, std::string> > sources_;
    std::unique_ptr<VirtualPolymeraseReader> currentReader_;

private:
    void OpenNextReader(void);
};

} // namespace BAM
} // namespace PacBio

#endif // VIRTUALPOLYMERASECOMPOSITEREADER_H
