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
// Author: Derek Barnett

#ifndef PULSE2BASECACHE_H
#define PULSE2BASECACHE_H

#include "pbbam/Config.h"
#include <boost/dynamic_bitset.hpp>
#include <string>
#include <cassert>
#include <cctype>
namespace PacBio {
namespace BAM {
namespace internal {

class Pulse2BaseCache
{
public:
    /// \brief Creates a Pulse2BaseCache from pulseCall data ('pc' tag)
    ///
    /// Computes & stores cache of basecalled vs. squashed pulse positions for
    /// later masking of pulse data.
    ///
    /// \param pulseCalls[in]   string contents of 'pc' tag
    ///
    Pulse2BaseCache(const std::string& pulseCalls)
        : data_(pulseCalls.size())
    {
        // basecalled pulse -> data[i] == 1
        // squashed pulse   -> data[i] == 0
        //
        const auto numPulses = pulseCalls.size();
        for (size_t i = 0; i < numPulses; ++i)
            data_[i] = std::isupper(pulseCalls.at(i));
    }

    Pulse2BaseCache(void) = delete;
    Pulse2BaseCache(const Pulse2BaseCache& other) = default;
    Pulse2BaseCache(Pulse2BaseCache&& other) = default;
    Pulse2BaseCache& operator=(const Pulse2BaseCache&) = default;
    Pulse2BaseCache& operator=(Pulse2BaseCache&&) = default;
    ~Pulse2BaseCache(void) noexcept {}

public:

    ///
    /// \brief FindFirst
    /// \return
    ///
    size_t FindFirst(void) const
    { return data_.find_first(); }

    ///
    /// \brief FindNext
    /// \param from
    /// \return
    ///
    size_t FindNext(size_t from) const
    { return data_.find_next(from); }

    ///
    /// \brief IsBasecallAt
    /// \param pos
    /// \return
    ///
    bool IsBasecallAt(const size_t pos) const
    { return data_[pos]; }

    /// \returns the total number of pulses (basecalled & squashed)
    ///
    size_t NumPulses(void) const
    {
        return data_.size();
    }

    /// \returns the total number of basecalled pulses
    ///
    size_t NumBases(void) const
    {
        return data_.count();
    }

    /// \brief Removes squashed pulse positions from input data.
    ///
    /// \param[in]  Contents of any per-pulse tag.
    /// \returns    Input \p pulseData less all squashed pulses
    ///
    template<typename T>
    T RemoveSquashedPulses(const T& pulseData) const
    {
        const auto numPulses = pulseData.size();
        assert(numPulses == data_.size());

        // The reserve() below overshoots the required space, but numPulses is cheap
        // to compute, and by definition will be sufficient to hold the result. Thus
        // we only ever need to do one allocation.
        //
        T result;
        result.reserve(numPulses);

        // Only include data at positions that match our cached pulse data.
        //
        size_t inputIndex = 0;
        for (size_t i = 0; i < numPulses; ++i) {
            if (data_[i])
                result.push_back(pulseData.at(inputIndex));
            ++inputIndex;
        }
        return result;
    }

private:
    boost::dynamic_bitset<> data_;
};

} // namespace internal
} // namespace BAM
} // namespace PacBio

#endif // PULSE2BASECACHE_H
