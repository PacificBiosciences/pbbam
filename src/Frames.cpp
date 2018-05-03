// File Description
/// \file Frames.cpp
/// \brief Implements the Frames class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/Frames.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <mutex>

namespace PacBio {
namespace BAM {
namespace internal {

static std::vector<uint16_t> framepoints;
static std::vector<uint8_t> frameToCode;
static uint16_t maxFramepoint;
static std::mutex initIpdDownsamplingMutex;

static void InitIpdDownsampling()
{
    std::lock_guard<std::mutex> lock(initIpdDownsamplingMutex);

    if (!framepoints.empty()) return;

    // liftover from Dave's python code:
    // .../bioinformatics/tools/kineticsTools/kineticsTools/_downsampling.py

    const int B = 2;
    const int t = 6;
    const double T = pow(B, t);

    int next = 0;
    double grain;
    const int end = 256 / T;
    for (int i = 0; i < end; ++i) {
        grain = pow(B, i);
        std::vector<uint16_t> nextOnes;
        for (double j = 0; j < T; ++j)
            nextOnes.push_back(j * grain + next);
        next = nextOnes.back() + grain;
        framepoints.insert(framepoints.end(), nextOnes.cbegin(), nextOnes.cend());
    }
    assert(framepoints.size() - 1 <= std::numeric_limits<uint8_t>::max());

    const uint16_t maxElement = (*max_element(framepoints.cbegin(), framepoints.cend()));
    frameToCode.assign(maxElement + 1, 0);

    const int fpEnd = framepoints.size() - 1;
    uint8_t i = 0;
    uint16_t fl = 0;
    uint16_t fu = 0;
    for (; i < fpEnd; ++i) {
        fl = framepoints[i];
        fu = framepoints[i + 1];
        if (fu > fl + 1) {
            const int middle = (fl + fu) / 2;
            for (int f = fl; f < middle; ++f)
                frameToCode[f] = i;
            for (int f = middle; f < fu; ++f)
                frameToCode[f] = i + 1;
        } else
            frameToCode[fl] = i;
    }

    // this next line differs from the python implementation (there, it's "i+1")
    // our C++ for loop has incremented our index counter one more time than the indexes from python enumerate(...)
    frameToCode[fu] = i;
    maxFramepoint = fu;
}

static inline uint16_t CodeToFrames(const uint8_t code) { return framepoints[code]; }

static std::vector<uint16_t> CodeToFrames(const std::vector<uint8_t>& codedData)
{
    InitIpdDownsampling();

    const auto length = codedData.size();
    std::vector<uint16_t> frames(length, 0);
    for (size_t i = 0; i < length; ++i)
        frames[i] = CodeToFrames(codedData[i]);
    return frames;
}

static inline uint8_t FramesToCode(const uint16_t frame)
{
    return frameToCode[std::min(maxFramepoint, frame)];
}

static std::vector<uint8_t> FramesToCode(const std::vector<uint16_t>& frames)
{
    InitIpdDownsampling();

    const auto length = frames.size();
    std::vector<uint8_t> result(length, 0);
    for (size_t i = 0; i < length; ++i)
        result[i] = FramesToCode(frames[i]);
    return result;
}

}  // namespace internal

Frames::Frames() {}

Frames::Frames(std::vector<uint16_t> frames) : data_{std::move(frames)} {}

Frames Frames::Decode(const std::vector<uint8_t>& codedData)
{
    return Frames{internal::CodeToFrames(codedData)};
}

std::vector<uint8_t> Frames::Encode(const std::vector<uint16_t>& frames)
{
    return internal::FramesToCode(frames);
}

}  // namespace BAM
}  // namespace PacBio
