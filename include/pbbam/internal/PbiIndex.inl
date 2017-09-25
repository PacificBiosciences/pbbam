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
/// \file PbiIndex.inl
/// \brief Inline implementations for the PbiIndex class.
//
// Author: Derek Barnett

#include "pbbam/BamRecord.h"
#include "pbbam/PbiFile.h"
#include "pbbam/PbiIndex.h"
#include "pbbam/PbiRawData.h"

#include <algorithm>
#include <map>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>
#include <cassert>

namespace PacBio {
namespace BAM {
namespace internal {

// --------------------------
// Pbi Lookup Aggregate
// --------------------------

class PbiIndexPrivate
{
public:
    PbiIndexPrivate() = default;
    PbiIndexPrivate(const PbiRawData& rawIndex);
    PbiIndexPrivate(PbiRawData&& rawIndex);

    PbiIndexPrivate(const PbiIndexPrivate&) = delete;
    PbiIndexPrivate& operator=(const PbiIndexPrivate&) = delete;

    std::unique_ptr<PbiIndexPrivate> DeepCopy() const;

public:
    bool HasSection(const PbiFile::Section flag) const;
    void SetSection(const PbiFile::Section flag, bool ok = true);

public:
    IndexResultBlocks LookupReference(const int32_t tId) const;

private:
    IndexResultBlocks MergeBlocksWithOffsets(const IndexList& indices) const;

public:
    std::string filename_;
    PbiFile::VersionEnum version_ = PbiFile::CurrentVersion;
    PbiFile::Sections sections_ = PbiFile::BASIC;
    uint32_t numReads_ = 0;

    // lookup structures
    BasicLookupData     basicData_;
    MappedLookupData    mappedData_;
    ReferenceLookupData referenceData_;
    BarcodeLookupData   barcodeData_;
};

inline bool PbiIndexPrivate::HasSection(const PbiFile::Section flag) const
{ return (sections_ & flag) != 0; }

inline void PbiIndexPrivate::SetSection(const PbiFile::Section flag, bool ok)
{ if (ok) sections_ |= flag; else sections_ &= ~flag; }

inline IndexResultBlocks
PbiIndexPrivate::LookupReference(const int32_t tId) const
{
    if (!HasSection(PbiFile::REFERENCE))
        return IndexResultBlocks{ };

    const auto& indexRange = referenceData_.Indices(tId);
    if (indexRange.first == nullIndex() && indexRange.second == nullIndex())
        return IndexResultBlocks{ };
    const auto numReads = indexRange.second - indexRange.first;
    auto blocks = IndexResultBlocks{ IndexResultBlock(indexRange.first, numReads) };
    basicData_.ApplyOffsets(blocks);
    return blocks;
}

inline IndexResultBlocks
PbiIndexPrivate::MergeBlocksWithOffsets(const IndexList& indices) const
{
    auto blocks = mergedIndexBlocks(indices);
    basicData_.ApplyOffsets(blocks);
    return blocks;
}

} // namespace internal

inline PbiFile::Sections PbiIndex::FileSections() const
{ return d_->sections_; }

inline bool PbiIndex::HasBarcodeData() const
{ return d_->HasSection(PbiFile::BARCODE); }

inline bool PbiIndex::HasMappedData() const
{ return d_->HasSection(PbiFile::MAPPED); }

inline bool PbiIndex::HasReferenceData() const
{ return d_->HasSection(PbiFile::REFERENCE); }

inline bool PbiIndex::HasSection(const PbiFile::Section section) const
{ return d_->HasSection(section); }

inline uint32_t PbiIndex::NumReads() const
{ return d_->numReads_; }

inline PbiFile::VersionEnum PbiIndex::Version() const
{ return d_->version_; }

inline const BarcodeLookupData& PbiIndex::BarcodeData() const
{ return d_->barcodeData_; }

inline const BasicLookupData& PbiIndex::BasicData() const
{ return d_->basicData_; }

inline const MappedLookupData& PbiIndex::MappedData() const
{ return d_->mappedData_; }

inline const ReferenceLookupData& PbiIndex::ReferenceData() const
{ return d_->referenceData_; }

} // namespace BAM
} // namespace PacBio
