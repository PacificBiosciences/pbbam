#ifndef PBBAM_PBIRAWDATA_H
#define PBBAM_PBIRAWDATA_H

#include <pbbam/Config.h>

#include <pbbam/PbiFile.h>

#include <memory>
#include <string>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace BAM {

class BamFile;
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

    /// \brief Creates an empty data structure, preallocating space for a known
    ///        number of records.
    PbiRawBarcodeData(std::uint32_t numReads);

    PbiRawBarcodeData() = default;

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

    std::vector<std::int16_t> bcForward_;
    std::vector<std::int16_t> bcReverse_;
    std::vector<std::int8_t> bcQual_;

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

    /// \brief Creates an empty data structure, preallocating space for a known
    ///        number of records.
    PbiRawMappedData(std::uint32_t numReads);

    PbiRawMappedData() = default;

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
    std::uint32_t NumDeletedBasesAt(std::size_t recordIndex) const;

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
    std::uint32_t NumInsertedBasesAt(std::size_t recordIndex) const;

    /// \brief Calculates the number of deleted & inserted bases for a
    ///        particular record.
    ///
    /// \param[in] recordIndex  i-th record in the data set
    /// \returns a pair consisting of (numDeletions,numInsertions)
    ///
    std::pair<std::uint32_t, std::uint32_t> NumDeletedAndInsertedBasesAt(
        std::size_t recordIndex) const;

    /// \}

public:
    /// \name Raw Data Containers
    /// \{

    std::vector<std::int32_t> tId_;
    std::vector<std::uint32_t> tStart_;
    std::vector<std::uint32_t> tEnd_;
    std::vector<std::uint32_t> aStart_;
    std::vector<std::uint32_t> aEnd_;
    std::vector<std::uint8_t> revStrand_;
    std::vector<std::uint32_t> nM_;
    std::vector<std::uint32_t> nMM_;
    std::vector<std::uint8_t> mapQV_;
    std::vector<std::uint32_t> nInsOps_;
    std::vector<std::uint32_t> nDelOps_;

    bool hasIndelOps_ = true;

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
    using ID = std::uint32_t;
    using Row = std::uint32_t;

public:
    static const ID UNMAPPED_ID;
    static const Row UNSET_ROW;

public:
    /// \name Constructors & Related Methods
    /// \{

    /// \brief Creates a default entry.
    ///
    /// - default ID:   PbiReferenceEntry::UNMAPPED_ID \n
    /// - default rows: PbiReferenceEntry::UNSET_ROW
    ///
    PbiReferenceEntry();

    /// \brief Creates a reference entry, with no rows set.
    ///
    /// - default rows: PbiReferenceEntry::UNSET_ROW
    ///
    PbiReferenceEntry(ID id);

    /// \brief Creates a reference entry, with rows set.
    ///
    PbiReferenceEntry(ID id, Row beginRow, Row endRow);

    bool operator==(const PbiReferenceEntry& other) const noexcept;

    /// \}

public:
    /// \name Reference Data Members
    /// \{

    ID tId_;
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

    /// \brief Creates an empty data structure, preallocating space for a
    ///        number of references.
    ///
    /// This constructor is recommended as this is the safest way to ensure that
    /// references without observed mappings are included in the final output.
    ///
    PbiRawReferenceData(std::uint32_t numRefs);

    PbiRawReferenceData() = default;

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

    /// \brief Creates an empty data structure, preallocating space for a known
    ///        number of records.
    PbiRawBasicData(std::uint32_t numReads);

    PbiRawBasicData() = default;

    /// \}

public:
    /// \name Index Construction
    /// \{

    /// \brief Adds a record's mapping data.
    ///
    /// \param[in] b        %BAM record
    /// \param[in] offset   \b virtual file offset where record begins
    ///
    void AddRecord(const BamRecord& b, std::int64_t offset);

    /// \}

public:
    /// \name Raw Data Containers
    /// \{

    std::vector<std::int32_t> rgId_;
    std::vector<std::int32_t> qStart_;
    std::vector<std::int32_t> qEnd_;
    std::vector<std::int32_t> holeNumber_;
    std::vector<float> readQual_;
    std::vector<std::uint8_t> ctxtFlag_;
    std::vector<std::int64_t> fileOffset_;
    std::vector<std::uint16_t> fileNumber_;

    /// \}
};

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

    /// \brief Loads raw PBI data from a file.
    ///
    /// \param[in] pbiFilename      ".pbi" filename
    ///
    /// \throws std::runtime_error if file contents cannot be loaded properly
    ///
    PbiRawData(std::string pbiFilename);

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

    PbiRawData() = default;

    /// \}

public:
    /// \name PBI General Attributes
    /// \{

    /// \returns true if index has BarcodeData section
    bool HasBarcodeData() const;

    /// \returns true if index has MappedData section
    bool HasMappedData() const;

    /// \returns true if index has ReferenceData section
    bool HasReferenceData() const;

    /// \returns true if index has \b section
    /// \param[in] section PbiFile::Section identifier
    ///
    bool HasSection(PbiFile::Section section) const;

    /// \returns index filename ("*.pbi")
    ///
    /// \note Returns an empty string if the underlying data was calculated in
    ///       code or aggregated from a DataSet, rather than loaded from a
    ///       single PBI file.
    ///
    std::string Filename() const;

    /// \returns enum flags representing the file sections present
    PbiFile::Sections FileSections() const;

    /// \returns the number of records in the PBI(s)
    std::uint32_t NumReads() const;

    /// \returns the PBI file's version
    PbiFile::VersionEnum Version() const;

    /// \}

public:
    /// \name Raw Data Components
    /// \{

    /// \returns const reference to BarcodeData lookup structure
    ///
    /// May be empty, check result of HasBarcodeData.
    ///
    const PbiRawBarcodeData& BarcodeData() const;

    /// \returns const reference to BasicData lookup structure
    const PbiRawBasicData& BasicData() const;

    /// \returns const reference to MappedData lookup structure
    ///
    /// May be empty, check result of HasMappedData.
    ///
    const PbiRawMappedData& MappedData() const;

    /// \returns const reference to reference data lookup structure
    ///
    /// May be empty, check result of HasReferenceData.
    ///
    const PbiRawReferenceData& ReferenceData() const;

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
    PbiRawData& NumReads(std::uint32_t num);

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
    PbiRawBarcodeData& BarcodeData();

    /// \returns reference to BasicData lookup structure
    PbiRawBasicData& BasicData();

    /// \returns reference to MappedData lookup structure
    ///
    /// May be empty, check result of HasMappedData.
    ///
    PbiRawMappedData& MappedData();

    /// \returns reference to reference data lookup structure
    ///
    /// May be empty, check result of HasReferenceData.
    ///
    PbiRawReferenceData& ReferenceData();

    /// \}

private:
    std::string filename_;
    PbiFile::VersionEnum version_ = PbiFile::CurrentVersion;
    PbiFile::Sections sections_ = PbiFile::ALL;
    std::uint32_t numReads_ = 0;
    PbiRawBarcodeData barcodeData_;
    PbiRawMappedData mappedData_;
    PbiRawReferenceData referenceData_;
    PbiRawBasicData basicData_;
};

// PBI index caching

using PbiIndexCache = std::shared_ptr<std::vector<std::shared_ptr<PbiRawData>>>;

PbiIndexCache MakePbiIndexCache(const DataSet& dataset);
PbiIndexCache MakePbiIndexCache(const std::vector<BamFile>&);
PbiIndexCache MakePbiIndexCache(const BamFile& bamFile);

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_PBIRAWDATA_H
