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
/// \file BamRecordImpl.inl
/// \brief Inline implementations for the BamRecordImpl class.
//
// Author: Derek Barnett

#include "pbbam/BamRecordImpl.h"

namespace PacBio {
namespace BAM {

inline uint32_t BamRecordImpl::Bin() const
{ return d_->core.bin; }

inline BamRecordImpl& BamRecordImpl::Bin(uint32_t bin)
{ d_->core.bin = bin; return *this; }

inline uint32_t BamRecordImpl::Flag() const
{ return d_->core.flag; }

inline BamRecordImpl& BamRecordImpl::Flag(uint32_t flag)
{ d_->core.flag = flag; return *this; }

inline int32_t BamRecordImpl::InsertSize() const
{ return d_->core.isize; }

inline BamRecordImpl& BamRecordImpl::InsertSize(int32_t iSize)
{ d_->core.isize = iSize; return *this; }

inline uint8_t BamRecordImpl::MapQuality() const
{ return d_->core.qual; }

inline BamRecordImpl& BamRecordImpl::MapQuality(uint8_t mapQual)
{ d_->core.qual = mapQual; return *this; }

inline PacBio::BAM::Position BamRecordImpl::MatePosition() const
{ return d_->core.mpos; }

inline BamRecordImpl& BamRecordImpl::MatePosition(PacBio::BAM::Position pos)
{ d_->core.mpos = pos; return *this; }

inline int32_t BamRecordImpl::MateReferenceId() const
{ return d_->core.mtid; }

inline BamRecordImpl& BamRecordImpl::MateReferenceId(int32_t id)
{ d_->core.mtid = id; return *this; }

inline PacBio::BAM::Position BamRecordImpl::Position() const
{ return d_->core.pos; }

inline BamRecordImpl& BamRecordImpl::Position(PacBio::BAM::Position pos)
{ d_->core.pos = pos; return *this; }

inline int32_t BamRecordImpl::ReferenceId() const
{ return d_->core.tid; }

inline BamRecordImpl& BamRecordImpl::ReferenceId(int32_t id)
{ d_->core.tid = id; return *this; }

inline bool BamRecordImpl::IsDuplicate() const
{ return (d_->core.flag & BamRecordImpl::DUPLICATE) != 0; }

inline BamRecordImpl& BamRecordImpl::SetDuplicate(bool ok)
{
    if (ok) d_->core.flag |=  BamRecordImpl::DUPLICATE;
    else    d_->core.flag &= ~BamRecordImpl::DUPLICATE;
    return *this;
}

inline bool BamRecordImpl::IsFailedQC() const
{ return (d_->core.flag & BamRecordImpl::FAILED_QC) != 0; }

inline BamRecordImpl& BamRecordImpl::SetFailedQC(bool ok)
{
    if (ok) d_->core.flag |=  BamRecordImpl::FAILED_QC;
    else    d_->core.flag &= ~BamRecordImpl::FAILED_QC;
    return *this;
}

inline bool BamRecordImpl::IsFirstMate() const
{ return (d_->core.flag & BamRecordImpl::MATE_1) != 0; }

inline BamRecordImpl& BamRecordImpl::SetFirstMate(bool ok)
{
    if (ok) d_->core.flag |=  BamRecordImpl::MATE_1;
    else    d_->core.flag &= ~BamRecordImpl::MATE_1;
    return *this;
}

inline bool BamRecordImpl::IsMapped() const
{ return (d_->core.flag & BamRecordImpl::UNMAPPED) == 0; }

inline BamRecordImpl& BamRecordImpl::SetMapped(bool ok)
{
    if (ok) d_->core.flag &= ~BamRecordImpl::UNMAPPED;
    else    d_->core.flag |=  BamRecordImpl::UNMAPPED;
    return *this;
}

inline bool BamRecordImpl::IsMateMapped() const
{ return (d_->core.flag & BamRecordImpl::MATE_UNMAPPED) == 0; }

inline BamRecordImpl& BamRecordImpl::SetMateMapped(bool ok)
{
    if (ok) d_->core.flag &= ~BamRecordImpl::MATE_UNMAPPED;
    else    d_->core.flag |=  BamRecordImpl::MATE_UNMAPPED;
    return *this;
}

inline bool BamRecordImpl::IsMateReverseStrand() const
{ return (d_->core.flag & BamRecordImpl::MATE_REVERSE_STRAND) != 0; }

inline BamRecordImpl& BamRecordImpl::SetMateReverseStrand(bool ok)
{
    if (ok) d_->core.flag |=  BamRecordImpl::MATE_REVERSE_STRAND;
    else    d_->core.flag &= ~BamRecordImpl::MATE_REVERSE_STRAND;
    return *this;
}

inline bool BamRecordImpl::IsPaired() const
{ return (d_->core.flag & BamRecordImpl::PAIRED) != 0; }

inline BamRecordImpl& BamRecordImpl::SetPaired(bool ok)
{
    if (ok) d_->core.flag |=  BamRecordImpl::PAIRED;
    else    d_->core.flag &= ~BamRecordImpl::PAIRED;
    return *this;
}

inline bool BamRecordImpl::IsPrimaryAlignment() const
{ return (d_->core.flag & BamRecordImpl::SECONDARY) == 0; }

inline BamRecordImpl& BamRecordImpl::SetPrimaryAlignment(bool ok)
{
    if (ok) d_->core.flag &= ~BamRecordImpl::SECONDARY;
    else    d_->core.flag |=  BamRecordImpl::SECONDARY;
    return *this;
}

inline bool BamRecordImpl::IsProperPair() const
{ return (d_->core.flag & BamRecordImpl::PROPER_PAIR) != 0; }

inline BamRecordImpl& BamRecordImpl::SetProperPair(bool ok)
{
    if (ok) d_->core.flag |=  BamRecordImpl::PROPER_PAIR;
    else    d_->core.flag &= ~BamRecordImpl::PROPER_PAIR;
    return *this;
}

inline bool BamRecordImpl::IsReverseStrand() const
{ return (d_->core.flag & BamRecordImpl::REVERSE_STRAND) != 0; }

inline BamRecordImpl& BamRecordImpl::SetReverseStrand(bool ok)
{
    if (ok) d_->core.flag |=  BamRecordImpl::REVERSE_STRAND;
    else    d_->core.flag &= ~BamRecordImpl::REVERSE_STRAND;
    return *this;
}

inline bool BamRecordImpl::IsSecondMate() const
{ return (d_->core.flag & BamRecordImpl::MATE_2) != 0; }

inline BamRecordImpl& BamRecordImpl::SetSecondMate(bool ok)
{
    if (ok) d_->core.flag |=  BamRecordImpl::MATE_2;
    else    d_->core.flag &= ~BamRecordImpl::MATE_2;
    return *this;
}

inline bool BamRecordImpl::IsSupplementaryAlignment() const
{ return (d_->core.flag & BamRecordImpl::SUPPLEMENTARY) != 0; }

inline BamRecordImpl& BamRecordImpl::SetSupplementaryAlignment(bool ok)
{
    if (ok) d_->core.flag |=  BamRecordImpl::SUPPLEMENTARY;
    else    d_->core.flag &= ~BamRecordImpl::SUPPLEMENTARY;
    return *this;
}

} // namespace BAM
} // namespace PacBio
