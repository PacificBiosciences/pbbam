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
/// \file BamRecordImpl.h
/// \brief Defines the BamRecordImpl class.
//
// Author: Derek Barnett

#ifndef BAMRECORDIMPL_H
#define BAMRECORDIMPL_H

#include "pbbam/BamRecordTag.h"
#include "pbbam/Cigar.h"
#include "pbbam/Config.h"
#include "pbbam/Position.h"
#include "pbbam/QualityValues.h"
#include "pbbam/TagCollection.h"
#include <htslib/sam.h>
#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <string>

namespace PacBio {
namespace BAM {

namespace internal { class BamRecordMemory; }

/// \brief The BamRecordImpl class holds all data necessary for creating,
///        querying or editing a generic %BAM record.
///
/// For PacBio-specific extensions and convenience methods, see BamRecord.
///
/// \note This class is mostly an internal implementation detail and will
///       likely be removed from the public API in the future. Please use
///       BamRecord as much as possible.
///
class PBBAM_EXPORT BamRecordImpl
{
public:

    /// These flags describe the alignment status of the record.
    enum AlignmentFlag
    {
        PAIRED              = 0x0001    ///< Record comes from paired-end sequencing
      , PROPER_PAIR         = 0x0002    ///< Each mate of a pair was properly aligned ("proper" as determined by aligner)
      , UNMAPPED            = 0x0004    ///< Record was not mapped by aligner
      , MATE_UNMAPPED       = 0x0008    ///< Record's mate was not mapped by aligner
      , REVERSE_STRAND      = 0x0010    ///< Record was aligned to reverse strand (Sequence() is reverse-complemented)
      , MATE_REVERSE_STRAND = 0x0020    ///< Record's mate was aligned to reverse strand (mate's Sequence() is reverse-complemented)
      , MATE_1              = 0x0040    ///< Record is first mate of pair
      , MATE_2              = 0x0080    ///< Record is second mate of pair
      , SECONDARY           = 0x0100    ///< Record is a secondary alignment
      , FAILED_QC           = 0x0200    ///< Record failed quality controls
      , DUPLICATE           = 0x0400    ///< Record is a PCR/optical duplicate
      , SUPPLEMENTARY       = 0x0800    ///< Record is a supplementary alignment
    };

public:
    /// \name Constructors & Related Methods
    /// \{

    BamRecordImpl();
    BamRecordImpl(const BamRecordImpl& other);
    BamRecordImpl(BamRecordImpl&& other);
    BamRecordImpl& operator=(const BamRecordImpl& other);
    BamRecordImpl& operator=(BamRecordImpl&& other);
    virtual ~BamRecordImpl() = default;

    /// \}

public:
    /// \name Core Data
    /// \{

    /// \returns this record's assigned (BAI) index bin ID.
    uint32_t Bin() const;

    /// \returns this record's alignment flag, in raw integer form.
    uint32_t Flag() const;

    /// \returns this record's insert size
    int32_t InsertSize() const;

    /// \returns this record's mapping quality. A value of 255 indicates "unknown"
    uint8_t MapQuality() const;

    /// \returns this record's mate's mapped position, or -1 if unmapped
    PacBio::BAM::Position MatePosition() const;

    /// \returns this record's mate's mapped reference ID, or -1 if unmapped
    int32_t MateReferenceId() const;

    /// \returns this record's mapped position, or -1 if unmapped
    PacBio::BAM::Position Position() const;

    /// \returns this record's mate's mapped reference ID, or -1 if unmapped
    int32_t ReferenceId() const;

    /// Sets the record's (BAI) index bin ID.
    ///
    /// \param[in] bin BAI index bin ID.
    /// \returns reference to this record
    ///
    BamRecordImpl& Bin(uint32_t bin);

    /// Sets this record's alignment flag, using a raw integer.
    ///
    /// \param[in] flag raw alignment flag
    /// \returns reference to this record
    ///
    BamRecordImpl& Flag(uint32_t flag);

    /// Sets this record's insert size.
    ///
    /// \param[in] iSize insert size
    /// \returns reference to this record
    ///
    BamRecordImpl& InsertSize(int32_t iSize);

    /// Sets this record's map quality.
    ///
    /// \param[in] mapQual mapping quality - value of 255 indicates "unknown"
    /// \returns reference to this record
    ///
    BamRecordImpl& MapQuality(uint8_t mapQual);

    /// Sets this record's mate's mapped position.
    ///
    /// \param[in] pos mapped position. A value of -1 indicates unmapped.
    /// \returns reference to this record
    ///
    BamRecordImpl& MatePosition(PacBio::BAM::Position pos);

    /// Sets this record's mate's mapped reference ID
    ///
    /// \param[in] id reference ID. A value of -1 indicates unmapped.
    /// \returns reference to this record
    ///
    BamRecordImpl& MateReferenceId(int32_t id);

    /// Sets this record's mapped position.
    ///
    /// \param[in] pos mapped position. A value of -1 indicates unmapped.
    /// \returns reference to this record
    ///
    BamRecordImpl& Position(PacBio::BAM::Position pos);

    /// Sets this record's mapped reference ID
    ///
    /// \param[in] id reference ID. A value of -1 indicates unmapped.
    /// \returns reference to this record
    ///
    BamRecordImpl& ReferenceId(int32_t id);

    /// \}

public:
    /// \name Alignment Flags
    /// \{

    /// \returns true if this record is a PCR/optical duplicate
    bool IsDuplicate() const;

    /// \returns true if this record failed quality controls
    bool IsFailedQC() const;

    /// \returns true if this record is the first mate of a pair
    bool IsFirstMate() const;

    /// \returns true if this record was mapped by aligner
    bool IsMapped() const;

    /// \returns true if this record's mate was mapped by aligner
    bool IsMateMapped() const;

    /// \returns true if this record's mate was mapped to the reverse strand
    bool IsMateReverseStrand() const;

    /// \returns true if this record comes from paired-end sequencing
    bool IsPaired() const;

    /// \returns true if this record is a read's primary alignment
    bool IsPrimaryAlignment() const;

    /// \returns true if this record & its mate were properly aligned
    bool IsProperPair() const;

    /// \returns true if this record was mapped to the reverse strand
    bool IsReverseStrand() const;

    /// \returns true if this record is the second mate of a pair
    bool IsSecondMate() const;

    /// \returns true if this record is a supplementary alignment
    bool IsSupplementaryAlignment() const;

    /// Sets whether this record is a PCR/optical duplicate
    BamRecordImpl& SetDuplicate(bool ok);

    /// Sets whether this record failed quality controls
    BamRecordImpl& SetFailedQC(bool ok);

    /// Sets whether this record is the first mate of a pair.
    BamRecordImpl& SetFirstMate(bool ok);

    /// Sets whether this record was aligned.
    BamRecordImpl& SetMapped(bool ok);

    /// Sets whether this record's mate was aligned.
    BamRecordImpl& SetMateMapped(bool ok);

    /// Sets whether this record's mate mapped to reverse strand.
    BamRecordImpl& SetMateReverseStrand(bool ok);

    /// Sets whether this record came from paired-end sequencing.
    BamRecordImpl& SetPaired(bool ok);

    /// Sets whether this record is a read's primary alignment.
    BamRecordImpl& SetPrimaryAlignment(bool ok);

    /// Sets whether this record & its mate were properly mapped, per the aligner.
    BamRecordImpl& SetProperPair(bool ok);

    /// Sets whether this record mapped to reverse strand.
    BamRecordImpl& SetReverseStrand(bool ok);

    /// Sets whether this record is the second mate of a pair.
    BamRecordImpl& SetSecondMate(bool ok);

    /// Sets whether this record is a supplementary alignment.
    BamRecordImpl& SetSupplementaryAlignment(bool ok);

    /// \}

public:
    /// \name Variable-length Data (sequence, qualities, etc.)
    /// \{

    /// \returns the record's CIGAR data as a Cigar object
    Cigar CigarData() const;

    /// Sets the record's CIGAR data using a Cigar object
    ///
    /// \param[in] cigar PacBio::BAM::Cigar object
    /// \returns reference to this record
    ///
    BamRecordImpl& CigarData(const Cigar& cigar);

    /// Sets the record's CIGAR data using a CIGAR-formatted string.
    ///
    /// \param[in] cigarString CIGAR-formatted string
    /// \returns reference to this record
    ///
    BamRecordImpl& CigarData(const std::string& cigarString);

    // TODO: CIGAR iterator - Cigar only or here as well ??

    /// \returns the record's query name
    std::string Name() const;

    /// Sets the record's "query name".
    ///
    /// \param name new name
    /// \returns reference to this record
    ///
    BamRecordImpl& Name(const std::string& name);

    /// \returns the record's quality values (phred-style ASCII)
    ///
    /// \note Usually Qualities().size() == Sequence.size(). However, in
    ///       some data sets, the quality values are not provided. In that
    ///       case, this method will return an empty container.
    ///
    QualityValues Qualities() const;

    /// \returns the record's DNA sequence.
    std::string Sequence() const;

    size_t SequenceLength() const;

    /// \brief Sets the record's DNA sequence and quality values
    ///
    /// This is an overloaded function. Sets the DNA sequence and quality
    /// values, using the length of \p sequence.
    ///
    /// \note When using this overload (and \p qualities is non-empty), the
    ///       lengths of \p sequence and \p qualities \b must be equal.
    ///
    /// \todo How to handle mismatched lengths?
    ///
    /// \param[in] sequence  std::string containing DNA sequence
    /// \param[in] qualities std::string containing ASCII quality values
    ///
    /// \returns reference to this record.
    ///
    /// \sa SetSequenceAndQualities(const char* sequence,
    ///     const size_t sequenceLength, const char* qualities)
    ///
    BamRecordImpl& SetSequenceAndQualities(const std::string& sequence,
                                           const std::string& qualities = std::string());

    /// \brief Sets the record's DNA sequence and quality values.
    ///
    /// The \p sequence must consist of IUPAC nucleotide codes {=ACMGRSVTWYHKDBN}.
    /// The \p qualities, if not empty, must consist of 'phred'-style ASCII
    /// quality values. \p qualities may be an empty string or NULL pointer in
    /// cases where there are no such data available.
    ///
    /// \param[in] sequence         C-string containing DNA sequence
    /// \param[in] sequenceLength   length of DNA sequence
    /// \param[in] qualities        C-string containing 'phred-style' ASCII
    ///                             quality values
    ///
    /// \note \p sequence does \b NOT have to be NULL-terminated. Length is
    ///       explicitly determined by the value of \p sequenceLength provided.
    ///
    /// \returns reference to this record.
    ///
    BamRecordImpl& SetSequenceAndQualities(const char* sequence,
                                           const size_t sequenceLength,
                                           const char* qualities = nullptr);

    /// \brief Sets the record's DNA sequence and quality values.
    ///
    /// The \p encodedSequence should be preencoded/packed into the BAM binary
    /// format. The \p qualities, if not empty, must consist of 'phred'-style
    /// ASCII quality values. \p qualities may be an empty string or NULL
    /// pointer in cases where there are no such data available.
    ///
    /// \param[in] encodedSequence      C-string containing BAM-format-encoded
    ///                                 DNA sequence
    /// \param[in] rawSequenceLength    length of DNA sequence (not the encoded
    ///                                 length)
    /// \param[in] qualities            C-string containing 'phred-style' ASCII
    ///                                 quality values
    ///
    /// \note \p encodedSequence does \b NOT have to be NULL-terminated. Length
    ///       is explicitly determined by the value of \p sequenceLength
    ///       provided.
    ///
    /// \returns reference to this record.
    ///
    /// \sa SetSequenceAndQualities(const char* sequence,
    ///     const size_t sequenceLength, const char* qualities)
    ///
    BamRecordImpl& SetPreencodedSequenceAndQualities(const char* encodedSequence,
                                                     const size_t rawSequenceLength,
                                                     const char* qualities = nullptr);

    /// \}

public:
    /// \name Tag Data
    /// \{

    /// \returns record's full tag data as a TagCollection object
    TagCollection Tags() const;

    /// \brief Sets the record's full tag data via a TagCollection object
    ///
    BamRecordImpl& Tags(const TagCollection& tags);

    /// \brief Adds a new tag to this record.
    ///
    /// \param[in] tagName  2-character tag name.
    /// \param[in] value    Tag object that describes the type & value of data
    ///                     to be added
    ///
    /// \note Any value that can be used to implicitly construct a Tag is valid.
    /// \code
    ///     string s;
    ///     vector<uint32_t> v;
    ///     record.AddTag("XX", s); // will add a string-type tag
    ///     record.AddTag("YY", v); // will add a uint32-array-type tag
    /// \endcode
    ///
    /// \returns true if tag was successfully added.
    ///
    bool AddTag(const std::string& tagName,
                const Tag& value);

    /// \brief Adds a new tag to this record.
    ///
    /// This is an overloaded method.
    ///
    /// \param[in] tag      BamRecordTag enum
    /// \param[in] value    Tag object that describes the type & value of data
    ///                     to be added
    /// \returns true if tag was successfully added.
    ///
    bool AddTag(const BamRecordTag tag,
                const Tag& value);

    /// \brief Adds a new tag to this record, with an optional modifier.
    ///
    /// \param[in] tagName              2-character tag name.
    /// \param[in] value                Tag object that describes the type &
    ///                                 value of data to be added
    /// \param[in] additionalModifier   optional extra modifier (for explicit
    ///                                 modification of an otherwise const Tag)
    ///
    /// \note Any value that can be used to implicitly construct a Tag is valid.
    /// \code
    ///     char c;
    ///     string h;
    ///     record.AddTag("XX", c, TagModifier::ASCII_CHAR); // will add a char-type tag
    ///     record.AddTag("YY", h, TagModifier::HEX_STRING); // will add a hex string-type tag
    /// \endcode
    ///
    /// \returns true if tag was successfully added.
    ///
    bool AddTag(const std::string& tagName,
                const Tag& value,
                const TagModifier additionalModifier);

    /// \brief Adds a new tag to this record, with an optional modifier.
    ///
    /// This is an overloaded method.
    ///
    /// \param[in] tag                  BamRecordTag enum.
    /// \param[in] value                Tag object that describes the type &
    ///                                 value of data to be added
    /// \param[in] additionalModifier   optional extra modifier (for explicit
    ///                                 modification of an otherwise const Tag)
    ///
    /// \returns true if tag was successfully added.
    ///
    bool AddTag(const BamRecordTag tag,
                const Tag& value,
                const TagModifier additionalModifier);

    /// \brief Edits an existing tag on this record.
    ///
    /// \param[in] tagName      2-character tag name. Name must be present
    ///                         (see HasTag)
    /// \param[in] newValue     Tag object that describes the type & value of
    ///                         new data to be added
    ///
    /// \note Any value that can be used to implicitly construct a Tag is valid.
    /// \code
    ///     string s;
    ///     vector<uint32_t> v;
    ///     record.EditTag("XX", s); // will overwrite tag XX with a string-type tag
    ///     record.EditTag("YY", v); // will overwrite tag YY with a uint32-array-type tag
    /// \endcode
    ///
    /// \returns true if tag was successfully edited.
    ///
    bool EditTag(const std::string& tagName,
                 const Tag& newValue);

    /// \brief Edits an existing tag on this record.
    ///
    /// This is an overloaded method.
    ///
    /// \param[in] tag          BamRecordTag enum
    /// \param[in] newValue     Tag object that describes the type & value of
    ///                         new data to be added
    ///
    /// \returns true if tag was successfully edited.
    ///
    bool EditTag(const BamRecordTag tag,
                 const Tag& newValue);

    /// \brief Edits an existing tag on this record.
    ///
    /// \param[in] tagName              2-character tag name. Name must be
    ///                                 present (see HasTag)
    /// \param[in] value                Tag object that describes the type &
    ///                                 value of new data to be added
    /// \param[in] additionalModifier   optional extra modifier (for explicit
    ///                                 modification of an otherwise const Tag)
    ///
    /// \note Any value that can be used to implicitly construct a Tag is valid.
    /// \code
    ///     char c;
    ///     string h;
    ///     record.EditTag("XX", c, TagModifier::ASCII_CHAR); // will overwrite tag XX with a char-type tag
    ///     record.EditTag("YY", h, TagModifier::HEX_STRING); // will overwrite tag YY with a hex string-type tag
    /// \endcode
    ///
    /// \returns true if tag was successfully edited.
    ///
    bool EditTag(const std::string& tagName,
                 const Tag& value,
                 const TagModifier additionalModifier);

    /// \brief Edits an existing tag on this record.
    ///
    /// This is an overloaded method.
    ///
    /// \param[in] tag                  BamRecordTag enum
    /// \param[in] value                Tag object that describes the type &
    ///                                 value of new data to be added
    /// \param[in] additionalModifier   optional extra modifier (for explicit
    ///                                 modification of an otherwise const Tag)
    ///
    /// \returns true if tag was successfully edited.
    ///
    bool EditTag(const BamRecordTag tag,
                 const Tag& value,
                 const TagModifier additionalModifier);


    /// \returns true if a tag with this name is present in this record.
    bool HasTag(const std::string& tagName) const;

    /// \returns true if this tag is present in this record.
    ///
    /// This is an overloaded method.
    ///
    bool HasTag(const BamRecordTag tag) const;

    /// \brief Removes an existing tag from this record.
    ///
    /// \param[in] tagName  2-character tag name.
    ///
    /// \returns true if tag was actaully removed (i.e. false if tagName
    ///          previously unknown)
    /// \sa HasTag
    ///
    bool RemoveTag(const std::string& tagName);

    /// \brief Removes an existing tag from this record.
    ///
    /// This is an overloaded method.
    ///
    /// \param[in] tag  BamRecordTag enum
    ///
    /// \returns true if tag was actaully removed (i.e. false if tagName
    ///          previously unknown)
    /// \sa HasTag
    ///
    bool RemoveTag(const BamRecordTag tag);

    /// \brief Fetches a tag from this record.
    ///
    /// \param[in] tagName  2-character tag name.
    ///
    /// \returns Tag object for the requested name. If name is unknown, a
    ///          default constructed Tag is returned (Tag::IsNull() is true).
    ///
    Tag TagValue(const std::string& tagName) const;

    /// \brief Fetches a tag from this record.
    ///
    /// This is an overloaded method
    ///
    /// \param[in] tag  BamRecordTag enum
    ///
    /// \returns Tag object for the requested name. If name is unknown, a
    ///          default constructed Tag is returned (Tag::IsNull() is true).
    ///
    Tag TagValue(const BamRecordTag tag) const;

    // change above to Tag();

//    template<typename T>
//    T TagValue(const std::string& tagName) const;


    /// \}

private:
    // returns a BamRecordImpl object, with a deep copy of @rawData contents
    static BamRecordImpl FromRawData(const std::shared_ptr<bam1_t>& rawData);

    // internal memory setup/expand methods
    void InitializeData();
    void MaybeReallocData();
    void UpdateTagMap() const; // allowed to be called from const methods
                                   // (lazy update on request)

    // internal tag helper methods
    bool AddTagImpl(const std::string& tagName,
                    const Tag& value,
                    const TagModifier additionalModifier);
    bool RemoveTagImpl(const std::string& tagName);
    int TagOffset(const std::string& tagName) const;

    // core seq/qual logic shared by the public API
    BamRecordImpl& SetSequenceAndQualitiesInternal(const char* sequence,
                                                      const size_t sequenceLength,
                                                      const char* qualities,
                                                      bool isPreencoded);

private:

    // data members
    std::shared_ptr<bam1_t> d_;
    mutable std::map<uint16_t, int> tagOffsets_;

    // friends
    friend class internal::BamRecordMemory;
};

} // namespace BAM
} // namespace PacBio

#include "pbbam/internal/BamRecordImpl.inl"

#endif // BAMRECORDIMPL_H
