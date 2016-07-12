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
/// \file PbiRawData.h
/// \brief Defines the classes used for working with raw PBI data.
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
class DataSet;

/// \brief The PbiRawBarcodeData class represents the raw data stored in the
///        "BarcodeData" section of the PBI index.
///
class PBBAM_EXPORT PbiRawBarcodeData
{
public:
    /// \name Constructors & Related Methods
    /// \{

    /// \brief Creates an empty data structure.
    PbiRawBarcodeData(void);

    /// \brief Creates an empty data structure, preallocating space for a known
    ///        number of records.
    PbiRawBarcodeData(uint32_t numReads);

    PbiRawBarcodeData(const PbiRawBarcodeData& other);
    PbiRawBarcodeData(PbiRawBarcodeData&& other);
    PbiRawBarcodeData& operator=(const PbiRawBarcodeData& other);
    PbiRawBarcodeData& operator=(PbiRawBarcodeData&& other);

    /// \}

public:
    /// \name Index Construction
    /// \{

    /// \brief Adds a record's barcode data.
    ///
    /// \param[in] b    %BAM record
    ///
    void AddRecord(const BamRecord& b);

    /// \}

public:
    /// \name Raw Data Containers
    /// \{

    std::vector<int16_t> bcForward_;
    std::vector<int16_t> bcReverse_;
    std::vector<int8_t>  bcQual_;

    /// \}
};

/// \brief The PbiRawMappedData class represents the raw data stored in the
///        "MappedData" section of the PBI index.
///
class PBBAM_EXPORT PbiRawMappedData
{
public:
    /// \name Constructors & Related Methods
    /// \{

    /// \brief Creates an empty data structure.
    PbiRawMappedData(void);

    /// \brief Creates an empty data structure, preallocating space for a known
    ///        number of records.
    PbiRawMappedData(uint32_t numReads);

    PbiRawMappedData(const PbiRawMappedData& other);
    PbiRawMappedData(PbiRawMappedData&& other);
    PbiRawMappedData& operator=(const PbiRawMappedData& other);
    PbiRawMappedData& operator=(PbiRawMappedData&& other);

    /// \}

public:
    /// \name Index Construction
    /// \{

    /// \brief Adds a record's mapping data.
    ///
    /// \param[in] b    %BAM record
    ///
    void AddRecord(const BamRecord& b);

    /// \}

public:
    /// \name Index Data Query
    /// \{

    /// \brief Calculates the number of deleted bases for a particular record.
    ///
    /// Convenvience method. Equivalent to:
    /// \code{.cpp}
    /// NumDeletedAndInsertedBasesAt(i).first;
    /// \endcode
    ///
    /// \param[in] recordIndex  i-th record
    /// \returns number of deleted bases
    ///
    uint32_t NumDeletedBasesAt(size_t recordIndex) const;

    /// \brief Calculates the number of inserted bases for a particular record.
    ///
    /// Convenvience method. Equivalent to:
    /// \code{.cpp}
    /// NumDeletedAndInsertedBasesAt(i).second;
    /// \endcode
    ///
    /// \param[in] recordIndex  i-th record
    /// \returns number of inserted bases
    ///
    uint32_t NumInsertedBasesAt(size_t recordIndex) const;

    /// \brief Calculates the number of deleted & inserted bases for a
    ///        particular record.
    ///
    /// \param[in] recordIndex  i-th record in the data set
    /// \returns a pair consisting of (numDeletions,numInsertions)
    ///
    std::pair<uint32_t, uint32_t>
    NumDeletedAndInsertedBasesAt(size_t recordIndex) const;

    /// \}

public:
    /// \name Raw Data Containers
    /// \{

    std::vector<int32_t>  tId_;
    std::vector<uint32_t> tStart_;
    std::vector<uint32_t> tEnd_;
    std::vector<uint32_t> aStart_;
    std::vector<uint32_t> aEnd_;
    std::vector<uint8_t>  revStrand_;
    std::vector<uint32_t> nM_;
    std::vector<uint32_t> nMM_;
    std::vector<uint8_t>  mapQV_;

    /// \}
};

/// \brief The PbiReferenceEntryClass represents a single reference in the PBI
///        CoordinateSorted section.
///
/// A reference entry consists of an associated reference ID (tId), as well as
/// start and end indices into the %BAM or PBI.
///
/// \note Rows are given in the interval [start, end).
///
class PBBAM_EXPORT PbiReferenceEntry
{
public:
    typedef uint32_t ID;
    typedef uint32_t Row;

public:
    static const ID  UNMAPPED_ID;
    static const Row UNSET_ROW;

public:
    /// \name Constructors & Related Methods
    /// \{

    /// \brief Creates a default entry.
    ///
    /// - default ID:   PbiReferenceEntry::UNMAPPED_ID \n
    /// - default rows: PbiReferenceEntry::UNSET_ROW
    ///
    PbiReferenceEntry(void);

    /// \brief Creates a reference entry, with no rows set.
    ///
    /// - default rows: PbiReferenceEntry::UNSET_ROW
    ///
    PbiReferenceEntry(ID id);

    /// \brief Creates a reference entry, with rows set.
    ///
    PbiReferenceEntry(ID id, Row beginRow, Row endRow);

    PbiReferenceEntry(const PbiReferenceEntry& other);
    PbiReferenceEntry(PbiReferenceEntry&& other);
    PbiReferenceEntry& operator=(const PbiReferenceEntry& other);
    PbiReferenceEntry& operator=(PbiReferenceEntry&& other);

    bool operator==(const PbiReferenceEntry& other) const;

    /// \}

public:
    /// \name Reference Data Members
    /// \{

    ID  tId_;
    Row beginRow_;
    Row endRow_;

    /// \}
};

/// \brief The PbiRawReferenceData class represents the raw data stored in the
///        "CoordinateSortedData" section of the PBI index.
///
class PBBAM_EXPORT PbiRawReferenceData
{
public:
    /// \name Constructors & Related Methods
    /// \{

    /// \brief Creates an empty data structure.
    PbiRawReferenceData(void);

    /// \brief Creates an empty data structure, preallocating space for a
    ///        number of references.
    ///
    /// This constructor is recommended as this is the safest way to ensure that
    /// references without observed mappings are included in the final output.
    ///
    PbiRawReferenceData(uint32_t numRefs);

    PbiRawReferenceData(const PbiRawReferenceData& other);
    PbiRawReferenceData(PbiRawReferenceData&& other);
    PbiRawReferenceData& operator=(const PbiRawReferenceData& other);
    PbiRawReferenceData& operator=(PbiRawReferenceData&& other);

    /// \}

public:
    /// \name Raw Data Containers
    /// \{

    std::vector<PbiReferenceEntry> entries_;

    /// \}
};

/// \brief The PbiRawBasicData class represents the raw data stored in the
///        "BasicData" section of the PBI index.
///
class PBBAM_EXPORT PbiRawBasicData
{
public:
    /// \name Constructors & Related Methods
    /// \{

    /// \brief Creates an empty data structure.
    PbiRawBasicData(void);

    /// \brief Creates an empty data structure, preallocating space for a known
    ///        number of records.
    PbiRawBasicData(uint32_t numReads);

    PbiRawBasicData(const PbiRawBasicData& other);
    PbiRawBasicData(PbiRawBasicData&& other);
    PbiRawBasicData& operator=(const PbiRawBasicData& other);
    PbiRawBasicData& operator=(PbiRawBasicData&& other);

    /// \}

public:
    /// \name Index Construction
    /// \{

    /// \brief Adds a record's mapping data.
    ///
    /// \param[in] b        %BAM record
    /// \param[in] offset   \b virtual file offset where record begins
    ///
    void AddRecord(const BamRecord& b, int64_t offset);

    /// \}

public:
    /// \name Raw Data Containers
    /// \{

    std::vector<int32_t>  rgId_;
    std::vector<int32_t>  qStart_;
    std::vector<int32_t>  qEnd_;
    std::vector<int32_t>  holeNumber_;
    std::vector<float>    readQual_;
    std::vector<uint8_t>  ctxtFlag_;
    std::vector<int64_t>  fileOffset_;
    std::vector<uint16_t> fileNumber_;

    /// \}
};

/// \deprecated For legacy-code support only, and will be removed soon.
///             Use PbiRawBasicData instead.
///
typedef PbiRawBasicData PbiRawSubreadData;

/// \brief The PbiRawData class provides an representation of raw PBI index
///        data, used mostly for construction or I/O.
///
/// The PbiRawData class itself provides access to a few high-level attributes
/// (e.g. version, number of records, etc.). The actual index data is stored
/// in its member components:
///     PbiRawBasicData,
///     PbiRawMappedData,
///     PbiRawReferenceData, &
///     PbiRawBarcodeData .
///
class PBBAM_EXPORT PbiRawData
{
public:
    /// \name Constructors & Related Methods
    /// \{

    /// \brief Creates an empty raw data structure, ready for building.
    ///
    PbiRawData(void);

    /// \brief Loads raw PBI data from a file.
    ///
    /// \param[in] pbiFilename      ".pbi" filename
    ///
    /// \throws std::runtime_error if file contents cannot be loaded properly
    ///
    PbiRawData(const std::string& pbiFilename);

    /// \brief Loads a raw, aggregate PBI data from a dataset
    ///
    /// This constructor creates a raw index object that contains an aggregation
    /// of index data across the dataset.
    ///
    /// \note ReferenceData (the per-reference table for coordinate-sorted data)
    ///       is not currently available for the index aggregate. All other
    ///       per-record data sections will be present.
    ///
    /// \param[in] dataset  DataSet object
    ///
    /// \throws std::runtime_error if file(s) contents cannot be loaded properly
    ///
    explicit PbiRawData(const DataSet& dataset);

    PbiRawData(const PbiRawData& other);
    PbiRawData(PbiRawData&& other);
    PbiRawData& operator=(const PbiRawData& other);
    PbiRawData& operator=(PbiRawData&& other);
    ~PbiRawData(void);

    /// \}

public:
    /// \name PBI General Attributes
    /// \{

    /// \returns true if index has BarcodeData section
    bool HasBarcodeData(void) const;

    /// \returns true if index has MappedData section
    bool HasMappedData(void) const;

    /// \returns true if index has ReferenceData section
    bool HasReferenceData(void) const;

    /// \returns true if index has \b section
    /// \param[in] section PbiFile::Section identifier
    ///
    bool HasSection(const PbiFile::Section section) const;

    /// \returns index filename ("*.pbi")
    ///
    /// \note Returns an empty string if the underlying data was calculated in
    ///       code or aggregated from a DataSet, rather than loaded from a
    ///       single PBI file.
    ///
    std::string Filename(void) const;

    /// \returns enum flags representing the file sections present
    PbiFile::Sections FileSections(void) const;

    /// \returns the number of records in the PBI(s)
    uint32_t NumReads(void) const;

    /// \returns the PBI file's version
    PbiFile::VersionEnum Version(void) const;

    /// \}

public:
    /// \name Raw Data Components
    /// \{

    /// \returns const reference to BarcodeData lookup structure
    ///
    /// May be empty, check result of HasBarcodeData.
    ///
    const PbiRawBarcodeData& BarcodeData(void) const;

    /// \returns const reference to BasicData lookup structure
    const PbiRawBasicData& BasicData(void) const;

    /// \returns const reference to MappedData lookup structure
    ///
    /// May be empty, check result of HasMappedData.
    ///
    const PbiRawMappedData& MappedData(void) const;

    /// \returns const reference to reference data lookup structure
    ///
    /// May be empty, check result of HasReferenceData.
    ///
    const PbiRawReferenceData& ReferenceData(void) const;

    /// \}

public:
    /// \name PBI General Attributes
    /// \{

    /// \brief Sets the file section flags.
    ///
    /// \param[in] sections     section flags
    /// \returns reference to this index
    ///
    PbiRawData& FileSections(PbiFile::Sections sections);

    /// \brief Sets the number of indexed records.
    ///
    /// \param[in] num  number of records
    /// \returns reference to this index
    ///
    PbiRawData& NumReads(uint32_t num);

    /// \brief Sets PBI file version.
    ///
    /// \param[in] version  file version
    /// \returns reference to this index
    ///
    PbiRawData& Version(PbiFile::VersionEnum version);

    /// \}

public:
    /// \name Raw Data Components
    /// \{

    /// \returns reference to BarcodeData lookup structure
    ///
    /// May be empty, check result of HasBarcodeData.
    ///
    PbiRawBarcodeData& BarcodeData(void);

    /// \returns reference to BasicData lookup structure
    PbiRawBasicData& BasicData(void);

    /// \returns reference to MappedData lookup structure
    ///
    /// May be empty, check result of HasMappedData.
    ///
    PbiRawMappedData& MappedData(void);

    /// \returns reference to reference data lookup structure
    ///
    /// May be empty, check result of HasReferenceData.
    ///
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

} // namespace BAM
} // namespace PacBio

#include "pbbam/internal/PbiRawData.inl"

#endif // PBIRAWDATA_H
