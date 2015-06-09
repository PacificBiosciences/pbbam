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

// Author: Armin Töpfer

#ifndef VIRTUALPOLYMERASEREADER_H
#define VIRTUALPOLYMERASEREADER_H

#include <memory>

#include "pbbam/BamFile.h"
#include "pbbam/BamRecord.h"
#include "pbbam/Config.h"
#include "pbbam/EntireFileQuery.h"
#include "pbbam/virtual/VirtualPolymeraseBamRecord.h"

namespace PacBio {
namespace BAM {

class VirtualPolymeraseReader
{
public:
	/// Constructor takes two input bam file paths.
	/// \param[in] primaryBamFilePath hqregion.bam or subreads.bam file path
	/// \param[in] scrapsBamFilePath scraps.bam file path
    VirtualPolymeraseReader(const std::string& primaryBamFilePath,
                            const std::string& scrapsBamFilePath);

    VirtualPolymeraseReader() = delete;
    // Move constructor
    VirtualPolymeraseReader(VirtualPolymeraseReader&& src) = delete;
    // Copy constructor
    VirtualPolymeraseReader(const VirtualPolymeraseReader& src) = delete;
    // Move assignment operator
    VirtualPolymeraseReader& operator=(VirtualPolymeraseReader&& rhs) = delete;
    // Copy assignment operator
    VirtualPolymeraseReader& operator=(const VirtualPolymeraseReader& rhs) = delete;
    // Destructor
    ~VirtualPolymeraseReader() = default;

public:
	/// Provides the next stitched polymerase read
	VirtualPolymeraseBamRecord Next();

	/// Provides the next set of reads that belong to one ZMW.
	/// Enables stitching records in a distinct thread.
	std::vector<BamRecord> NextRaw();

	/// Returns true if more ZMWs are available for reading.
	bool HasNext();

	BamHeader PrimaryHeader();
	BamHeader ScrapsHeader();

private:
	const std::string                primaryBamFilePath_;
	const std::string                scrapsBamFilePath_;

	std::unique_ptr<BamFile>         primaryBamFile_;
	std::unique_ptr<BamFile>         scrapsBamFile_;
	std::unique_ptr<EntireFileQuery> primaryQuery_;
	std::unique_ptr<EntireFileQuery> scrapsQuery_;

	QueryBase::iterator              primaryIt_;
	QueryBase::iterator              scrapsIt_;

	std::unique_ptr<BamHeader>       polyHeader_;
};

} // namespace BAM
} // namespace PacBio

#endif // VIRTUALPOLYMERASEREADER_H
