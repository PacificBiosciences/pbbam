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
/// \file AlignmentPrinter.cpp
/// \brief Implements the AlignmentPrinter class.
//
// Author: Armin Töpfer

#include "pbbam/AlignmentPrinter.h"

#include <cmath>
#include <iostream>
#include <iomanip>  
#include <stdexcept>
#include <sstream>

using namespace PacBio;
using namespace PacBio::BAM;

AlignmentPrinter::AlignmentPrinter(const IndexedFastaReader& ifr)
    : ifr_(std::unique_ptr<IndexedFastaReader>(new IndexedFastaReader(ifr)))
{ }

std::string AlignmentPrinter::Print(const BamRecord& record,
                                    const Orientation orientation)
{
	std::string seq = record.Sequence(orientation, true, true);
	std::string ref = ifr_->ReferenceSubsequence(record, orientation, true, true);

	if (seq.size() != ref.size())
        throw std::runtime_error("Sequence and reference parts are of different size");

	int seqLength = 0;
	float matches = 0;
	std::string pretty;
    Position refCoord = record.ReferenceStart();
    Position seqCoord = record.QueryStart();

	for (size_t i = 0; i < seq.size();)
	{
		std::string refCoordStr = std::to_string(refCoord);
		std::string seqCoordStr = std::to_string(seqCoord);

		size_t maxCoordLength = std::max(refCoordStr.size(), seqCoordStr.size());
		while (refCoordStr.size() < maxCoordLength)
			refCoordStr = " "+refCoordStr;
		while (seqCoordStr.size() < maxCoordLength)
			seqCoordStr = " "+seqCoordStr;

		std::string seqWrap = seqCoordStr + " : ";
		std::string refWrap = refCoordStr + " : ";
		std::string prettyWrap(maxCoordLength+3, ' ');
		prettyWrap.reserve(seq.size());
		for (int j = 0; i < seq.size() && j < 40; ++i, ++j)
		{
			refWrap +=  ref[i];

			if (seq[i] == ref[i])
			{
				++matches;
				if (refCoord == 0 || refCoord % 10)
					prettyWrap += '|'; 
				else
				{
					prettyWrap += "\033[1m\x1b[31m";
					prettyWrap += '|'; 
					prettyWrap += "\033[0m\x1b[39;49m";
				}
				seqWrap += seq[i];
			}
			else if (seq[i] == '-' || ref[i] == '-')
			{
				prettyWrap +=  ' ';
				seqWrap += seq[i];
			}
			else
			{
				prettyWrap +=  '.';
				seqWrap += "\033[1m\x1b[31m";
				seqWrap += seq[i];
				seqWrap += "\033[0m\x1b[39;49m";
			}
			if (seq[i] != '-')
			{
				++seqLength;
				++seqCoord;
			}
			if (ref[i] != '-')
			{
				++refCoord;
			}
		}

		refCoordStr = std::to_string(refCoord);
		seqCoordStr = std::to_string(seqCoord);

		maxCoordLength = std::max(refCoordStr.size(), seqCoordStr.size());
		while (refCoordStr.size() < maxCoordLength)
			refCoordStr = " "+refCoordStr;
		while (seqCoordStr.size() < maxCoordLength)
			seqCoordStr = " "+seqCoordStr;

		seqWrap += " : " + seqCoordStr;
		refWrap += " : " + refCoordStr;

		pretty += refWrap + '\n' + prettyWrap + '\n' + seqWrap + "\n\n";
	}
	float similarity = matches/seq.size();

	std::stringstream output;

	output << "Read        : " << record.FullName() << std::endl;
	output << "Reference   : " << record.ReferenceName() << std::endl;
	output << std::endl;
	output << "Read-length : " << seqLength << std::endl;
	output << "Concordance : " << std::setprecision(3) << (similarity);
	output << std::endl;
	output << std::endl;
	output << pretty;

    return output.str();
}
