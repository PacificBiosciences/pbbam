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
/// \file PbiRawData.inl
/// \brief Inline implementations for the classes used for working with raw PBI
///        data.
//
// Author: Derek Barnett

#include "pbbam/PbiRawData.h"

namespace PacBio {
namespace BAM {

inline const PbiRawBarcodeData& PbiRawData::BarcodeData() const
{ return barcodeData_; }

inline PbiRawBarcodeData& PbiRawData::BarcodeData()
{ return barcodeData_; }

inline const PbiRawBasicData& PbiRawData::BasicData() const
{ return basicData_; }

inline PbiRawBasicData& PbiRawData::BasicData()
{ return basicData_; }

inline std::string PbiRawData::Filename() const
{ return filename_; }

inline PbiFile::Sections PbiRawData::FileSections() const
{ return sections_; }

inline PbiRawData& PbiRawData::FileSections(PbiFile::Sections sections)
{ sections_ = sections; return *this; }

inline bool PbiRawData::HasBarcodeData() const
{ return HasSection(PbiFile::BARCODE); }

inline bool PbiRawData::HasMappedData() const
{ return HasSection(PbiFile::MAPPED); }

inline bool PbiRawData::HasReferenceData() const
{ return HasSection(PbiFile::REFERENCE); }

inline bool PbiRawData::HasSection(const PbiFile::Section section) const
{ return (sections_ & section) != 0; }

inline uint32_t PbiRawData::NumReads() const
{ return numReads_; }

inline PbiRawData& PbiRawData::NumReads(uint32_t num)
{ numReads_ = num; return *this; }

inline const PbiRawMappedData& PbiRawData::MappedData() const
{ return mappedData_; }

inline PbiRawMappedData& PbiRawData::MappedData()
{ return mappedData_; }

inline const PbiRawReferenceData& PbiRawData::ReferenceData() const
{ return referenceData_; }

inline PbiRawReferenceData& PbiRawData::ReferenceData()
{ return referenceData_; }

inline PbiFile::VersionEnum PbiRawData::Version() const
{ return version_; }

inline PbiRawData& PbiRawData::Version(PbiFile::VersionEnum version)
{ version_ = version; return *this; }

inline bool PbiReferenceEntry::operator==(const PbiReferenceEntry& other) const
{
    return tId_      == other.tId_ &&
           beginRow_ == other.beginRow_ &&
           endRow_   == other.endRow_;
}

} // namespace BAM
} // namespace PacBio
