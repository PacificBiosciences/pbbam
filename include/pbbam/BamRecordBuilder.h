// File Description
/// \file BamRecordBuilder.h
/// \brief Defines the BamRecordBuilder class.
//
// Author: Derek Barnett

#ifndef BAMRECORDBUILDER_H
#define BAMRECORDBUILDER_H

#include <cstdint>
#include <string>
#include "pbbam/BamHeader.h"
#include "pbbam/BamRecord.h"
#include "pbbam/Config.h"

namespace PacBio {
namespace BAM {

/// \brief The BamRecordBuilder class provides a helper utility for building
///        BamRecords.
///
/// This class provides a mechanism for building up %BAM data and
/// lazy-encoding/constructing the actual BamRecord. Currently, the methods here
/// really only support  filling in the low-level SAM/BAM-style fields, not so
/// much the PacBio-specific fields.
///
class PBBAM_EXPORT BamRecordBuilder
{
public:
    /// \name Constructors & Related Methods
    /// \{

    /// \brief Creates an empty %BAM record builder.
    BamRecordBuilder();

    /// \brief Creates an empty %BAM record builder, with header info to apply
    ///        to built records.
    ///
    /// \param[in] header   BamHeader object
    ///
    explicit BamRecordBuilder(BamHeader header);

    /// \brief Creates record builder with inital record data.
    ///
    /// \param[in] prototype    data from this record will be used to seed the
    ///                         builder
    ///
    BamRecordBuilder(const BamRecord& prototype);

    BamRecordBuilder(const BamRecordBuilder&) = default;
    BamRecordBuilder(BamRecordBuilder&&) = default;
    BamRecordBuilder& operator=(const BamRecordBuilder&) = default;
    BamRecordBuilder& operator=(BamRecordBuilder&&) = default;
    ~BamRecordBuilder() = default;

    /// \}

public:
    /// \name Record-Building
    /// \{

    /// \brief Builds a BamRecord from current builder attributes.
    ///
    /// \returns newly-built BamRecord object
    ///
    BamRecord Build() const;

    /// \brief Replaces an existing BamRecord's data with current builder
    ///        attributes.
    ///
    /// \param[out] record resulting record
    /// \returns true if successful
    ///
    bool BuildInPlace(BamRecord& record) const;

    /// \brief Resets builder attributes to default values.
    ///
    void Reset();

    /// \brief Resets builder attributes with \p prototype's data.
    ///
    /// \param[in] prototype
    ///
    void Reset(BamRecord prototype);

    /// \}

public:
    /// \name Core Attribute Setup
    /// \{

    /// \brief Sets the record's (BAI) index bin ID.
    ///
    /// \param[in] bin BAI index bin ID.
    /// \returns reference to this builder
    ///
    BamRecordBuilder& Bin(const uint32_t bin);

    /// \brief Sets this record's alignment flag, using a raw integer.
    ///
    /// \param[in] flag raw alignment flag
    /// \returns reference to this record
    ///
    BamRecordBuilder& Flag(const uint32_t flag);

    /// \brief Sets this record's insert size.
    ///
    /// \param[in] iSize insert size
    /// \returns reference to this record
    ///
    BamRecordBuilder& InsertSize(const int32_t iSize);

    /// \brief Sets this record's map quality.
    ///
    /// \param[in] mapQual mapping quality - value of 255 indicates "unknown"
    /// \returns reference to this record
    ///
    BamRecordBuilder& MapQuality(const uint8_t mapQual);

    /// \brief Sets this record's mate's mapped position.
    ///
    /// \param[in] pos mapped position. A value of -1 indicates unmapped.
    /// \returns reference to this record
    ///
    BamRecordBuilder& MatePosition(const int32_t pos);

    /// \brief Sets this record's mate's mapped reference ID
    ///
    /// \param[in] id reference ID. A value of -1 indicates unmapped.
    /// \returns reference to this record
    ///
    BamRecordBuilder& MateReferenceId(const int32_t id);

    /// \brief Sets this record's mapped position.
    ///
    /// \param[in] pos mapped position. A value of -1 indicates unmapped.
    /// \returns reference to this record
    ///
    BamRecordBuilder& Position(const int32_t pos);

    /// \brief Sets this record's mapped reference ID
    ///
    /// \param[in] id reference ID. A value of -1 indicates unmapped.
    /// \returns reference to this record
    ///
    BamRecordBuilder& ReferenceId(const int32_t id);

    /// \}

public:
    /// \name Alignment Flag Setup
    /// \{

    /// \brief Sets whether this record is a PCR/optical duplicate
    BamRecordBuilder& SetDuplicate(bool ok);

    /// \brief Sets whether this record failed quality controls
    BamRecordBuilder& SetFailedQC(bool ok);

    /// \brief Sets whether this record is the first mate of a pair.
    BamRecordBuilder& SetFirstMate(bool ok);

    /// \brief Sets whether this record was aligned.
    BamRecordBuilder& SetMapped(bool ok);

    /// \brief Sets whether this record's mate was aligned.
    BamRecordBuilder& SetMateMapped(bool ok);

    /// \brief Sets whether this record's mate mapped to reverse strand.
    BamRecordBuilder& SetMateReverseStrand(bool ok);

    /// \brief Sets whether this record came from paired-end sequencing.
    BamRecordBuilder& SetPaired(bool ok);

    /// \brief Sets whether this record is a read's primary alignment.
    BamRecordBuilder& SetPrimaryAlignment(bool ok);

    /// \brief Sets whether this record & its mate were properly mapped, per the
    ///        aligner.
    ///
    BamRecordBuilder& SetProperPair(bool ok);

    /// \brief Sets whether this record mapped to reverse strand.
    BamRecordBuilder& SetReverseStrand(bool ok);

    /// \brief Sets whether this record is the second mate of a pair.
    BamRecordBuilder& SetSecondMate(bool ok);

    /// \brief Sets whether this record is a supplementary alignment.
    BamRecordBuilder& SetSupplementaryAlignment(bool ok);

    /// \}

public:
    /// \name Variable-Length Data Setup
    /// \{

    /// \brief Sets the record's CIGAR data.
    ///
    /// \returns reference to this builder
    ///
    BamRecordBuilder& Cigar(PacBio::BAM::Cigar cigar);

    /// \brief Sets the record's name.
    ///
    /// \returns reference to this builder
    ///
    BamRecordBuilder& Name(std::string name);

    /// \brief Sets the record's qualities.
    ///
    /// \returns reference to this builder
    ///
    BamRecordBuilder& Qualities(std::string qualities);

    /// \brief Sets the record's sequence.
    ///
    /// \returns reference to this builder
    ///
    BamRecordBuilder& Sequence(std::string sequence);

    /// \brief Sets the record's tags.
    ///
    /// \returns reference to this builder
    ///
    BamRecordBuilder& Tags(TagCollection tags);

    /// \}

private:
    BamHeader header_;
    bam1_core_t core_;
    std::string name_;
    std::string sequence_;
    std::string qualities_;
    PacBio::BAM::Cigar cigar_;
    TagCollection tags_;
};

}  // namespace BAM
}  // namespace PacBio

#include "pbbam/internal/BamRecordBuilder.inl"

#endif  // BAMRECORDBUILDER_H
