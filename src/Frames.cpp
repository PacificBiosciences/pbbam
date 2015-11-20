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
/// \file Frames.cpp
/// \brief Implements the Frames class.
//
// Author: Derek Barnett

#include "pbbam/Frames.h"
#include <algorithm>
#include <cassert>
#include <cmath>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

static vector<uint16_t> framepoints;
static vector<uint8_t> frameToCode;
static uint16_t maxFramepoint;

static
void InitIpdDownsampling(void)
{
    if (!framepoints.empty())
        return;

    // liftover from Dave's python code:
    // .../bioinformatics/tools/kineticsTools/kineticsTools/_downsampling.py

    const int B = 2;
    const int t = 6;
    const double T = pow(B, t);

    int next = 0;
    double grain;
    const int end = 256/T;
    for (int i = 0; i < end; ++i) {
        grain = pow(B, i);
        vector<uint16_t> nextOnes;
        for (double j = 0; j < T; ++j)
            nextOnes.push_back(j*grain + next);
        next = nextOnes.back() + grain;
        framepoints.insert(framepoints.end(), nextOnes.cbegin(), nextOnes.cend());
    }
    assert(framepoints.size()-1 <= UINT8_MAX);

    const uint16_t maxElement = (*max_element(framepoints.cbegin(), framepoints.cend()));
    frameToCode.assign(maxElement+1, 0);

    const int fpEnd = framepoints.size() - 1;
    uint8_t i = 0;
    uint16_t fl = 0;
    uint16_t fu = 0;
    for (; i < fpEnd; ++i) {
        fl = framepoints[i];
        fu = framepoints[i+1];
        if (fu > fl+1) {
            const int middle = (fl+fu)/2;
            for (int f = fl; f < middle; ++f)
                frameToCode[f] = i;
            for (int f = middle; f < fu; ++f)
                frameToCode[f] = i+1;
        } else
            frameToCode[fl] = i;
    }

    // this next line differs from the python implementation (there, it's "i+1")
    // our C++ for loop has incremented our index counter one more time than the indexes from python enumerate(...)
    frameToCode[fu] = i;
    maxFramepoint = fu;
}

static inline
uint16_t CodeToFrames(const uint8_t code)
{
    return framepoints[code];
}

static
vector<uint16_t> CodeToFrames(const vector<uint8_t>& codedData)
{
    InitIpdDownsampling();

    const size_t length = codedData.size();
    vector<uint16_t> frames(length, 0);
    for (size_t i = 0; i < length; ++i)
        frames[i] = CodeToFrames(codedData[i]);
    return frames;
}

static inline
uint8_t FramesToCode(const uint16_t frame)
{
    return frameToCode[std::min(maxFramepoint, frame)];
}

static
vector<uint8_t> FramesToCode(const vector<uint16_t>& frames)
{
    InitIpdDownsampling();

    const size_t length = frames.size();
    vector<uint8_t> result(length, 0);
    for (size_t i = 0; i < length; ++i)
        result[i] = FramesToCode(frames[i]);
    return result;
}

} // namespace internal
} // namespace BAM
} // namespace PacBio

Frames::Frames(void)
{ }

Frames::Frames(const std::vector<uint16_t>& frames)
    : data_(frames)
{ }

Frames::Frames(std::vector<uint16_t>&& frames)
    : data_(std::move(frames))
{ }

Frames::Frames(const Frames& other)
    : data_(other.data_)
{ }

Frames::Frames(Frames&& other)
    : data_(std::move(other.data_))
{ }

Frames::~Frames(void) { }

Frames& Frames::operator=(const Frames& other)
{ data_ = other.data_; return *this; }

Frames& Frames::operator=(Frames&& other)
{ data_ = std::move(other.data_); return *this; }

Frames Frames::Decode(const std::vector<uint8_t>& codedData)
{ return Frames(std::move(internal::CodeToFrames(codedData))); }

std::vector<uint8_t> Frames::Encode(const std::vector<uint16_t>& frames)
{ return internal::FramesToCode(frames); }
