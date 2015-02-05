// Copyright (c) 2014, Pacific Biosciences of California, Inc.
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

#ifndef SAMREADGROUP_H
#define SAMREADGROUP_H

#include "pbbam/Config.h"
#include <string>

namespace PacBio {
namespace BAM {

struct PBBAM_EXPORT SamReadGroup
{
public:
    SamReadGroup(void) { }
    SamReadGroup(const std::string& ID)
        : id(ID)
    { }
//    SamReadGroup(const SamReadGroup& other) = default;
//    SamReadGroup(SamReadGroup&& other) = default;
    ~SamReadGroup(void) { }

public:
    // DictionaryBase compatibility
    inline std::string Key(void) const;
    inline SamReadGroup& Key(const std::string& key);

public:
    std::string id;                     // ID:<Id>               * Unique ID required for valid SAM header*
    std::string sequencingCenter;       // CN:<SequencingCenter>
    std::string description;            // DS:<Description>
    std::string date;                   // DT:<Date>             * (ISO 8601)
    std::string flowOrder;              // FO:<FlowOrder>
    std::string keySequence;            // KS:<KeySequence>
    std::string library;                // LB:<Library>
    std::string programs;               // PG:<Programs>
    std::string predictedInsertSize;    // PI:<PredictedInsertSize>
    std::string platform;               // PL:<Platform>
    std::string platformUnit;           // PU:<PlatformUnit>
    std::string sample;                 // SM:<Sample>
};

inline std::string SamReadGroup::Key(void) const
{ return id; }

inline SamReadGroup& SamReadGroup::Key(const std::string& key)
{ id = key; return *this; }

} // namespace BAM
} // namespace PacBio

#endif // SAMREADGROUP_H
