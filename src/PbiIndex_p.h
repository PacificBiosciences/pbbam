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

#ifndef PBIINDEX_P_H
#define PBIINDEX_P_H

#include "pbbam/BamRecord.h"
#include "pbbam/PbiFile.h"
#include "pbbam/PbiIndex.h"
#include "pbbam/PbiRawData.h"

#include <map>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

namespace PacBio {
namespace BAM {
namespace internal {



//
// -----------------------------------------

struct LookupBase
{
    // general helper functions, typedefs, etc.
    typedef std::vector<size_t>       IndexList;
    typedef std::pair<size_t, size_t> IndexRange;

    // merging index ranges for result block



    static const size_t NullIndex;
};


struct PerReadLookupBase : public LookupBase
{
    // helper functions, typedefs, etc. for per-read info

    template<typename T>
    std::map<T, IndexList> MakeLookupMap(const std::vector<T>& rawData) const;

    template<typename T>
    std::map<T, IndexList> MakeLookupMap(std::vector<T>&& rawData) const;
};

struct SubreadLookupData : public PerReadLookupBase
{
    // ctors
    SubreadLookupData(void);
    SubreadLookupData(const PbiRawSubreadData& rawData);
    SubreadLookupData(PbiRawSubreadData&& rawData);

    // map ordering doesn't make sense, optimize for direct lookup
    std::unordered_map<int32_t, IndexList> rgId_;

    // numeric comparisons make sense, keep key ordering preserved
    std::map<int32_t,  IndexList> qStart_;
    std::map<int32_t,  IndexList> qEnd_;
    std::map<int32_t,  IndexList> holeNumber_;
    std::map<uint16_t, IndexList> readQual_;

    // offsets
    std::vector<int64_t> fileOffset_;
};

struct MappedLookupData : public PerReadLookupBase
{
    // ctors
    MappedLookupData(void);
    MappedLookupData(const PbiRawMappedData& rawData);
    MappedLookupData(PbiRawMappedData&& rawData);

    // numeric comparisons make sense, keep key ordering preserved
    std::map<int32_t,  IndexList> tId_;
    std::map<uint32_t, IndexList> tStart_;
    std::map<uint32_t, IndexList> tEnd_;
    std::map<uint32_t, IndexList> aStart_;
    std::map<uint32_t, IndexList> aEnd_;
    std::map<uint32_t, IndexList> nM_;
    std::map<uint32_t, IndexList> nMM_;
    std::map<uint8_t,  IndexList> mapQV_;

    // no need for map overhead, just store indices
    IndexList reverseStrand_;
    IndexList forwardStrand_;
};

struct ReferenceLookupData : public LookupBase
{
    // ctors
    ReferenceLookupData(void);
    ReferenceLookupData(const PbiRawReferenceData& rawData);
    ReferenceLookupData(PbiRawReferenceData&& rawData);

    // references_[tId] = (begin, end) indices
    // into SubreadLookupData::fileOffset_
    std::unordered_map<int32_t, IndexRange> references_;
};

struct BarcodeLookupData : public PerReadLookupBase
{
    // ctors
    BarcodeLookupData(void);
    BarcodeLookupData(const PbiRawBarcodeData& rawData);
    BarcodeLookupData(PbiRawBarcodeData&& rawData);

    // numeric comparisons make sense, keep key ordering preserved
    std::map<uint16_t, IndexList> bcLeft_;
    std::map<uint16_t, IndexList> bcRight_;
    std::map<uint8_t,  IndexList> bcQual_;

    // see if this works, or if can use 'direct' query ---|
    std::map<uint8_t, IndexList> ctxtFlag_;  // <---------|
};

class PbiIndexPrivate
{
public:
    PbiIndexPrivate(void);
    PbiIndexPrivate(const PbiRawData& rawIndex);
    PbiIndexPrivate(PbiRawData&& rawIndex);

    std::unique_ptr<PbiIndexPrivate> DeepCopy(void) const;


public:
    bool HasSection(const PbiFile::Section flag) const;
    void SetSection(const PbiFile::Section flag, bool ok = true);

public:
    PbiFile::VersionEnum version_;
    PbiFile::Sections sections_;
    uint32_t numReads_;

    // lookup structures
    SubreadLookupData   subreadData_;
    MappedLookupData    mappedData_;
    ReferenceLookupData referenceData_;
    BarcodeLookupData   barcodeData_;

private:
    // not-implemented - ensure no copy
    PbiIndexPrivate(const PbiIndexPrivate& other);
    PbiIndexPrivate& operator=(const PbiIndexPrivate& other);
};

inline bool PbiIndexPrivate::HasSection(const PbiFile::Section flag) const
{ return (sections_ & flag) != 0; }

inline void PbiIndexPrivate::SetSection(const PbiFile::Section flag, bool ok)
{ if (ok) sections_ |= flag; else sections_ &= ~flag; }

template<typename T> inline
std::map<T, PerReadLookupBase::IndexList>
PerReadLookupBase::MakeLookupMap(const std::vector<T>& rawData) const
{
    std::map<T, IndexList> result;
    const size_t numElements = rawData.size();
    for (size_t i = 0; i < numElements; ++i)
        result[ rawData.at(i) ].push_back(i);
    return result;
}

template<typename T> inline
std::map<T, PerReadLookupBase::IndexList>
PerReadLookupBase::MakeLookupMap(std::vector<T>&& rawData) const
{
    std::map<T, IndexList> result;
    const size_t numElements = rawData.size();
    for (size_t i = 0; i < numElements; ++i)
        result[ rawData.at(i) ].push_back(i);
    return result;
}

} // namespace internal
} // namespace BAM
} // namespace PacBio

#endif // PACBIOINDEX_P_H
