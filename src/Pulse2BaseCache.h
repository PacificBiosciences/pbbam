// Author: Derek Barnett

#ifndef PULSE2BASECACHE_H
#define PULSE2BASECACHE_H

#include <cassert>
#include <cctype>
#include <cstddef>
#include <string>

#include <boost/dynamic_bitset.hpp>

#include "pbbam/Config.h"

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
    Pulse2BaseCache(const std::string& pulseCalls) : data_(pulseCalls.size())
    {
        // basecalled pulse -> data[i] == 1
        // squashed pulse   -> data[i] == 0
        //
        const auto numPulses = pulseCalls.size();
        for (size_t i = 0; i < numPulses; ++i)
            data_[i] = std::isupper(pulseCalls.at(i));
    }

    Pulse2BaseCache() = delete;
    Pulse2BaseCache(const Pulse2BaseCache&) = default;
    Pulse2BaseCache(Pulse2BaseCache&&) = default;
    Pulse2BaseCache& operator=(const Pulse2BaseCache&) = default;
    Pulse2BaseCache& operator=(Pulse2BaseCache&&) = default;
    ~Pulse2BaseCache() = default;

public:
    ///
    /// \brief FindFirst
    /// \return
    ///
    size_t FindFirst() const { return data_.find_first(); }

    ///
    /// \brief FindNext
    /// \param from
    /// \return
    ///
    size_t FindNext(size_t from) const { return data_.find_next(from); }

    ///
    /// \brief IsBasecallAt
    /// \param pos
    /// \return
    ///
    bool IsBasecallAt(const size_t pos) const { return data_[pos]; }

    /// \returns the total number of pulses (basecalled & squashed)
    ///
    size_t NumPulses() const { return data_.size(); }

    /// \returns the total number of basecalled pulses
    ///
    size_t NumBases() const { return data_.count(); }

    /// \brief Removes squashed pulse positions from input data.
    ///
    /// \param[in]  Contents of any per-pulse tag.
    /// \returns    Input \p pulseData less all squashed pulses
    ///
    template <typename T>
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
            if (data_[i]) result.push_back(pulseData.at(inputIndex));
            ++inputIndex;
        }
        return result;
    }

private:
    boost::dynamic_bitset<> data_;
};

}  // namespace internal
}  // namespace BAM
}  // namespace PacBio

#endif  // PULSE2BASECACHE_H
