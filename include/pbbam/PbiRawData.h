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

#ifndef PBIRAWDATA_H
#define PBIRAWDATA_H

#include "pbbam/Config.h"
#include "pbbam/PbiFile.h"
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

class BamRecord;

class PBBAM_EXPORT PbiRawBarcodeData
{
public:
    PbiRawBarcodeData(void);
    PbiRawBarcodeData(uint32_t numReads);
    PbiRawBarcodeData(const PbiRawBarcodeData& other);
    PbiRawBarcodeData(PbiRawBarcodeData&& other);
    PbiRawBarcodeData& operator=(const PbiRawBarcodeData& other);
    PbiRawBarcodeData& operator=(PbiRawBarcodeData&& other);

public:
    /// Maybe add barcode data for \p b, if available.
    /// \returns true if record had barcode data
    ///
    bool AddRecord(const BamRecord& b);

public:
    std::vector<uint16_t> bcForward_;
    std::vector<uint16_t> bcReverse_;
    std::vector<uint8_t>  bcQual_;

};

class PBBAM_EXPORT PbiRawMappedData
{
public:
    PbiRawMappedData(void);
    PbiRawMappedData(uint32_t numReads);
    PbiRawMappedData(const PbiRawMappedData& other);
    PbiRawMappedData(PbiRawMappedData&& other);
    PbiRawMappedData& operator=(const PbiRawMappedData& other);
    PbiRawMappedData& operator=(PbiRawMappedData&& other);

public:
    /// Maybe add mapping data for \p b, if available.
    /// \returns true if record had mapping data
    ///
    bool AddRecord(const BamRecord& b);

    uint32_t NumDeletedBasesAt(size_t recordIndex) const;
    uint32_t NumInsertedBasesAt(size_t recordIndex) const;
    std::pair<uint32_t, uint32_t> NumDeletedAndInsertedBasesAt(size_t recordIndex) const;

public:
    std::vector<int32_t>  tId_;
    std::vector<uint32_t> tStart_;
    std::vector<uint32_t> tEnd_;
    std::vector<uint32_t> aStart_;
    std::vector<uint32_t> aEnd_;
    std::vector<uint8_t>  revStrand_;
    std::vector<uint32_t> nM_;
    std::vector<uint32_t> nMM_;
    std::vector<uint8_t>  mapQV_;
};

class PBBAM_EXPORT PbiReferenceEntry
{
public:
    typedef uint32_t ID;
    typedef uint32_t Row;

public:
    PbiReferenceEntry(void);
    PbiReferenceEntry(ID id);
    PbiReferenceEntry(ID id, Row beginRow, Row endRow);
    PbiReferenceEntry(const PbiReferenceEntry& other);
    PbiReferenceEntry(PbiReferenceEntry&& other);
    PbiReferenceEntry& operator=(const PbiReferenceEntry& other);
    PbiReferenceEntry& operator=(PbiReferenceEntry&& other);

    bool operator==(const PbiReferenceEntry& other) const
    {
        return tId_      == other.tId_ &&
               beginRow_ == other.beginRow_ &&
               endRow_   == other.endRow_;
    }

public:
    static const ID  UNMAPPED_ID;
    static const Row UNSET_ROW;

public:
    ID  tId_;
    Row beginRow_;
    Row endRow_;
};

class PBBAM_EXPORT PbiRawReferenceData
{
public:
    PbiRawReferenceData(void);
    PbiRawReferenceData(uint32_t numRefs);
    PbiRawReferenceData(const PbiRawReferenceData& other);
    PbiRawReferenceData(PbiRawReferenceData&& other);
    PbiRawReferenceData& operator=(const PbiRawReferenceData& other);
    PbiRawReferenceData& operator=(PbiRawReferenceData&& other);

public:
    std::vector<PbiReferenceEntry> entries_;
};

class PBBAM_EXPORT PbiRawBasicData
{
public:
    PbiRawBasicData(void);
    PbiRawBasicData(uint32_t numReads);
    PbiRawBasicData(const PbiRawBasicData& other);
    PbiRawBasicData(PbiRawBasicData&& other);
    PbiRawBasicData& operator=(const PbiRawBasicData& other);
    PbiRawBasicData& operator=(PbiRawBasicData&& other);

public:
    void AddRecord(const BamRecord& b, int64_t offset);

public:
    std::vector<int32_t>  rgId_;
    std::vector<int32_t>  qStart_;
    std::vector<int32_t>  qEnd_;
    std::vector<int32_t>  holeNumber_;
    std::vector<float>    readQual_;
    std::vector<uint8_t>  ctxtFlag_;
    std::vector<int64_t>  fileOffset_;
};

typedef PbiRawBasicData PbiRawSubreadData; // for existing code, possible in use externally

class PBBAM_EXPORT PbiRawData
{
public:
    /// \name Constructors & Related Methods
    /// \{

    /// Default ctor. Used in index building
    PbiRawData(void);

    /// Load raw data from \p pbiFilename.
    ///
    /// \param[in] pbiFilename PBI filename
    ///
    /// \throws if file contents cannot be loaded properly
    ///
    PbiRawData(const std::string& pbiFilename);

    PbiRawData(const PbiRawData& other);
    PbiRawData(PbiRawData&& other);
    PbiRawData& operator=(const PbiRawData& other);
    PbiRawData& operator=(PbiRawData&& other);
    ~PbiRawData(void);

    /// \}

public:
    /// \name Attributes
    /// \{

    bool HasBarcodeData(void) const;
    bool HasMappedData(void) const;
    bool HasReferenceData(void) const;
    bool HasSection(const PbiFile::Section section) const;

    /// \returns PBI filename (".pbi"), will be empty string if user-generated
    std::string Filename(void) const;
    PbiFile::Sections FileSections(void) const;
    uint32_t NumReads(void) const;
    PbiFile::VersionEnum Version(void) const;

    /// \}

public:
    /// \name Indexed Sections
    /// \{

    const PbiRawBarcodeData&   BarcodeData(void) const;
    const PbiRawBasicData&     BasicData(void) const;
    const PbiRawMappedData&    MappedData(void) const;
    const PbiRawReferenceData& ReferenceData(void) const;

    /// \}

public:
    /// \name Attributes
    /// \{

    PbiRawData& FileSections(PbiFile::Sections sections);
    PbiRawData& NumReads(uint32_t num);
    PbiRawData& Version(PbiFile::VersionEnum version);

    /// \}

public:
    /// \name Indexed Sections

    PbiRawBarcodeData&   BarcodeData(void);
    PbiRawBasicData&     BasicData(void);
    PbiRawMappedData&    MappedData(void);
    PbiRawReferenceData& ReferenceData(void);

    /// \}

private:
    std::string          filename_;
    PbiFile::VersionEnum version_;
    PbiFile::Sections    sections_;
    uint32_t             numReads_;
    PbiRawBarcodeData    barcodeData_;
    PbiRawMappedData     mappedData_;
    PbiRawReferenceData  referenceData_;
    PbiRawBasicData      basicData_;
};

inline const PbiRawBarcodeData& PbiRawData::BarcodeData(void) const
{ return barcodeData_; }

inline PbiRawBarcodeData& PbiRawData::BarcodeData(void)
{ return barcodeData_; }

inline std::string PbiRawData::Filename(void) const
{ return filename_; }

inline PbiFile::Sections PbiRawData::FileSections(void) const
{ return sections_; }

inline PbiRawData& PbiRawData::FileSections(PbiFile::Sections sections)
{ sections_ = sections; return *this; }

inline bool PbiRawData::HasBarcodeData(void) const
{ return HasSection(PbiFile::BARCODE); }

inline bool PbiRawData::HasMappedData(void) const
{ return HasSection(PbiFile::MAPPED); }

inline bool PbiRawData::HasReferenceData(void) const
{ return HasSection(PbiFile::REFERENCE); }

inline bool PbiRawData::HasSection(const PbiFile::Section section) const
{ return (sections_ & section) != 0; }

inline uint32_t PbiRawData::NumReads(void) const
{ return numReads_; }

inline PbiRawData& PbiRawData::NumReads(uint32_t num)
{ numReads_ = num; return *this; }

inline const PbiRawMappedData& PbiRawData::MappedData(void) const
{ return mappedData_; }

inline PbiRawMappedData& PbiRawData::MappedData(void)
{ return mappedData_; }

inline const PbiRawReferenceData& PbiRawData::ReferenceData(void) const
{ return referenceData_; }

inline PbiRawReferenceData& PbiRawData::ReferenceData(void)
{ return referenceData_; }

inline const PbiRawBasicData& PbiRawData::BasicData(void) const
{ return basicData_; }

inline PbiRawBasicData& PbiRawData::BasicData(void)
{ return basicData_; }

inline PbiFile::VersionEnum PbiRawData::Version(void) const
{ return version_; }

inline PbiRawData& PbiRawData::Version(PbiFile::VersionEnum version)
{ version_ = version; return *this; }

} // namespace BAM
} // namespace PacBio

#endif // PBIRAWDATA_H
