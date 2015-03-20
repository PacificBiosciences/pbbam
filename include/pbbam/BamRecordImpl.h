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

// Author: Derek Barnett

#ifndef BAMRECORDIMPL_H
#define BAMRECORDIMPL_H

#include "htslib/sam.h"
#include "pbbam/Cigar.h"
#include "pbbam/Config.h"
#include "pbbam/Position.h"
#include "pbbam/QualityValues.h"
#include "pbbam/TagCollection.h"
#include <memory>
#include <string>

namespace PacBio {
namespace BAM {

namespace internal { class BamRecordMemory; }

/// \brief The BamRecordImpl class holds all data necessary for creating,
/// querying or editing a valid, generic BAM record.
///
/// For PacBio-specific extensions and convenience methods, see BamRecord.
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

    BamRecordImpl(void);
    BamRecordImpl(const BamRecordImpl& other);
    BamRecordImpl(BamRecordImpl&& other);
    BamRecordImpl& operator=(const BamRecordImpl& other);
    BamRecordImpl& operator=(BamRecordImpl&& other);
    virtual ~BamRecordImpl(void);

    /// \}

public:

    /** @name Core Data
     */
    ///@{

    /// \returns this record's assigned (BAI) index bin ID.
    inline uint32_t Bin(void) const;

    /// \returns this record's alignment flag, in raw integer form.
    inline uint32_t Flag(void) const;

    /// \returns this record's insert size
    inline int32_t InsertSize(void) const;

    /// \returns this record's mapping quality. A value of 255 indicates "unknown"
    inline uint8_t MapQuality(void) const;

    /// \returns this record's mate's mapped position, or -1 if unmapped
    inline PacBio::BAM::Position MatePosition(void) const;

    /// \returns this record's mate's mapped reference ID, or -1 if unmapped
    inline int32_t MateReferenceId(void) const;

    /// \returns this record's mapped position, or -1 if unmapped
    inline PacBio::BAM::Position Position(void) const;

    /// \returns this record's mate's mapped reference ID, or -1 if unmapped
    inline int32_t ReferenceId(void) const;

    /// Sets the record's (BAI) index bin ID.
    ///
    /// \param[in] bin BAI index bin ID.
    /// \returns reference to this record
    inline BamRecordImpl& Bin(uint32_t bin);

    /// Sets this record's alignment flag, using a raw integer.
    ///
    /// \param[in] flag raw alignment flag
    /// \returns reference to this record
    inline BamRecordImpl& Flag(uint32_t flag);

    /// Sets this record's insert size.
    ///
    /// \param[in] iSize insert size
    /// \returns reference to this record
    inline BamRecordImpl& InsertSize(int32_t iSize);

    /// Sets this record's map quality.
    ///
    /// \param[in] mapQual mapping quality - value of 255 indicates "unknown"
    /// \returns reference to this record
    inline BamRecordImpl& MapQuality(uint8_t mapQual);

    /// Sets this record's mate's mapped position.
    ///
    /// \param[in] pos mapped position. A value of -1 indicates unmapped.
    /// \returns reference to this record
    inline BamRecordImpl& MatePosition(PacBio::BAM::Position pos);

    /// Sets this record's mate's mapped reference ID
    ///
    /// \param[in] id reference ID. A value of -1 indicates unmapped.
    /// \returns reference to this record
    inline BamRecordImpl& MateReferenceId(int32_t id);

    /// Sets this record's mapped position.
    ///
    /// \param[in] pos mapped position. A value of -1 indicates unmapped.
    /// \returns reference to this record
    inline BamRecordImpl& Position(PacBio::BAM::Position pos);

    /// Sets this record's mapped reference ID
    ///
    /// \param[in] id reference ID. A value of -1 indicates unmapped.
    /// \returns reference to this record
    inline BamRecordImpl& ReferenceId(int32_t id);

    ///@}

public:
    /** @name Alignment Flags
     */
    ///@{

    /// \returns true if this record is a PCR/optical duplicate
    inline bool IsDuplicate(void) const;

    /// \returns true if this record failed quality controls
    inline bool IsFailedQC(void) const;

    /// \returns true if this record is the first mate of a pair
    inline bool IsFirstMate(void) const;

    /// \returns true if this record was mapped by aligner
    inline bool IsMapped(void) const;

    /// \returns true if this record's mate was mapped by aligner
    inline bool IsMateMapped(void) const;

    /// \returns true if this record's mate was mapped to the reverse strand
    inline bool IsMateReverseStrand(void) const;

    /// \returns true if this record comes from paired-end sequencing
    inline bool IsPaired(void) const;

    /// \returns true if this record is a read's primary alignment
    inline bool IsPrimaryAlignment(void) const;

    /// \returns true if this record & its mate were properly aligned
    inline bool IsProperPair(void) const;

    /// \returns true if this record was mapped to the reverse strand
    inline bool IsReverseStrand(void) const;

    /// \returns true if this record is the second mate of a pair
    inline bool IsSecondMate(void) const;

    /// \returns true if this record is a supplementary alignment
    inline bool IsSupplementaryAlignment(void) const;

    /// Sets whether this record is a PCR/optical duplicate
    inline BamRecordImpl& SetDuplicate(bool ok);

    /// Sets whether this record failed quality controls
    inline BamRecordImpl& SetFailedQC(bool ok);

    /// Sets whether this record is the first mate of a pair.
    inline BamRecordImpl& SetFirstMate(bool ok);

    /// Sets whether this record was aligned.
    inline BamRecordImpl& SetMapped(bool ok);

    /// Sets whether this record's mate was aligned.
    inline BamRecordImpl& SetMateMapped(bool ok);

    /// Sets whether this record's mate mapped to reverse strand.
    inline BamRecordImpl& SetMateReverseStrand(bool ok);

    /// Sets whether this record came from paired-end sequencing.
    inline BamRecordImpl& SetPaired(bool ok);

    /// Sets whether this record is a read's primary alignment.
    inline BamRecordImpl& SetPrimaryAlignment(bool ok);

    /// Sets whether this record & its mate were properly mapped, per the aligner.
    inline BamRecordImpl& SetProperPair(bool ok);

    /// Sets whether this record mapped to reverse strand.
    inline BamRecordImpl& SetReverseStrand(bool ok);

    /// Sets whether this record is the second mate of a pair.
    inline BamRecordImpl& SetSecondMate(bool ok);

    /// Sets whether this record is a supplementary alignment.
    inline BamRecordImpl& SetSupplementaryAlignment(bool ok);

    ///@}

public:
    /** @name Variable-length Data (sequence, qualities, etc.)
     */
    ///@{

    /// \returns the record's CIGAR data as a Cigar object
    Cigar CigarData(void) const;

    /// Sets the record's CIGAR data using a Cigar object
    ///
    /// \param[in] cigar PacBio::BAM::Cigar object
    ///
    /// \returns reference to this record
    BamRecordImpl& CigarData(const Cigar& cigar);

    /// Sets the record's CIGAR data using a CIGAR-formatted string.
    ///
    /// \param[in] cigarString CIGAR-formatted string
    ///
    /// \returns reference to this record
    BamRecordImpl& CigarData(const std::string& cigarString);

    // TODO: CIGAR iterator - Cigar only or here as well ??

    /// \returns the record's query name
    std::string Name(void) const;

    /// Sets the record's "query name".
    ///
    /// \param name new name
    ///
    /// \returns reference to this record
    BamRecordImpl& Name(const std::string& name);

    /// \returns the record's quality values (phred-style ASCII)
    ///
    /// \note Usually Qualities().size() == Sequence.size(). However, in
    ///       some data sets, the quality values are not provided. In that
    ///       case, this method will return an empty container.
    QualityValues Qualities(void) const;

    /// \returns the record's DNA sequence.
    std::string Sequence(void) const;

    /// \brief Sets the record's DNA sequence and quality values
    ///
    /// This is an overloaded function. Sets the DNA sequence and quality values,
    /// using the length of \p sequence.
    ///
    /// \note When using this overload (and \p qualities is non-empty), the lengths
    ///       of \p sequence and \p qualities \b must be equal.
    ///
    /// \todo How to handle mismatched lenths?
    ///
    /// \param[in] sequence  std::string containing DNA sequence
    /// \param[in] qualities std::string containing ASCII quality values
    ///
    /// \returns reference to this record.
    ///
    /// \sa SetSequenceAndQualities(const char* sequence, const size_t sequenceLength, const char* qualities)
    ///
    BamRecordImpl& SetSequenceAndQualities(const std::string& sequence,
                                           const std::string& qualities = std::string());

    /// \brief Sets the record's DNA sequence and quality values.
    ///
    /// The \p sequence must consist of IUPAC nucleotide codes {=ACMGRSVTWYHKDBN}.
    /// The \p qualities, if not empty, must consist of 'phred'-style ASCII quality
    /// values. \p qualities may be an empty string or NULL pointer in cases where
    /// there are no such data available.
    ///
    /// \param[in] sequence       C-string containing DNA sequence
    /// \param[in] sequenceLength length of DNA sequence
    /// \param[in] qualities      C-string containing 'phred-style' ASCII quality values
    ///
    /// \note \p sequence does \b NOT have to be NULL-terminated. Length is explicitly
    ///        determined by the value of \p sequenceLength provided.
    ///
    /// \returns reference to this record.
    ///
    BamRecordImpl& SetSequenceAndQualities(const char* sequence,
                                           const size_t sequenceLength,
                                           const char* qualities = 0);

    /// \brief Sets the record's DNA sequence and quality values.
    ///
    /// The \p encodedSequence should be preencoded/packed into the BAM binary format.
    /// The \p qualities, if not empty, must consist of 'phred'-style ASCII quality values.
    /// \p qualities may be an empty string or NULL pointer in cases where there are no
    /// such data available.
    ///
    /// \param[in] encodedSequence   C-string containing BAM-format-encoded DNA sequence
    /// \param[in] rawSequenceLength length of DNA sequence (not the encoded length)
    /// \param[in] qualities         C-string containing 'phred-style' ASCII quality values
    ///
    /// \note \p encodedSequence does \b NOT have to be NULL-terminated. Length is explicitly
    ///        determined by the value of \p sequenceLength provided.
    ///
    /// \returns reference to this record.
    ///
    /// \sa SetSequenceAndQualities(const char* sequence, const size_t sequenceLength, const char* qualities)
    ///
    BamRecordImpl& SetPreencodedSequenceAndQualities(const char* encodedSequence,
                                                     const size_t rawSequenceLength,
                                                     const char* qualities = 0);

public:
    /** @name Tag Data
     */
    ///@{

    /// \returns record's full tag data as a TagCollection object
    TagCollection Tags(void) const;

    /// Sets the record's full tag data via a TagCollection object
    BamRecordImpl& Tags(const TagCollection& tags);

    /// Adds a new tag to this record.
    ///
    /// \param[in] tagName 2-character tag name.
    /// \param[in] value Tag object that describes the type & value of data to be added
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
    bool AddTag(const std::string& tagName, const Tag& value);

    /// Edits an existing tag on this record.
    ///
    /// \param[in] tagName 2-character tag name. Name must be present (see HasTag)
    /// \param[in] newValue Tag object that describes the type & value of new data to be added
    ///
    /// \note Any value that can be used to implicitly construct a Tag is valid.
    /// \code
    ///     string s;
    ///     vector<uint32_t> v;
    ///     record.EditTag("XX", s); // will overwrite tag XX with a string-type Tag
    ///     record.EditTag("YY", v); // will overwrite tag YY with a uint32-array-type Tag
    /// \endcode
    ///
    /// \returns true if tag was successfully edited.
    bool EditTag(const std::string& tagName, const Tag& newValue);

    /// \returns true if a tag with this name is present in this record.
    bool HasTag(const std::string& tagName) const;

    /// Removes an existing tag from this record.
    ///
    /// \param[in] tagName 2-character tag name.
    ///
    /// \returns true if tag was actaully removed (i.e. false if tagName previously unknown)
    /// \sa HasTag
    bool RemoveTag(const std::string& tagName);

    /// Fetches a tag from this record.
    ///
    /// \param[in] tagName 2-character tag name.
    ///
    /// \returns Tag object for the requested name. If name is unknown, a default constructed
    ///          Tag is returned (Tag::IsNull() is true).
    Tag TagValue(const std::string& tagName) const;

    // tag iterator   ?
    // tag operator[] ?

    ///@}

private:
    // returns a BamRecordImpl object, with a deep copy of @rawData contents
    static BamRecordImpl FromRawData(const std::shared_ptr<bam1_t>& rawData);

    // internal memory setup/expand methods
    void InitializeData(void);
    void MaybeReallocData(void);

    // core seq/qual logic shared by the public API
    BamRecordImpl& SetSequenceAndQualitiesInternal(const char* sequence,
                                                      const size_t sequenceLength,
                                                      const char* qualities,
                                                      bool isPreencoded);

private:

    // data members
    std::shared_ptr<bam1_t> d_;

    // friends
    friend class internal::BamRecordMemory;
};

inline uint32_t BamRecordImpl::Bin(void) const
{ return d_->core.bin; }

inline BamRecordImpl& BamRecordImpl::Bin(uint32_t bin)
{ d_->core.bin = bin; return *this; }

inline uint32_t BamRecordImpl::Flag(void) const
{ return d_->core.flag; }

inline BamRecordImpl& BamRecordImpl::Flag(uint32_t flag)
{ d_->core.flag = flag; return *this; }

inline int32_t BamRecordImpl::InsertSize(void) const
{ return d_->core.isize; }

inline BamRecordImpl& BamRecordImpl::InsertSize(int32_t iSize)
{ d_->core.isize = iSize; return *this; }

inline uint8_t BamRecordImpl::MapQuality(void) const
{ return d_->core.qual; }

inline BamRecordImpl& BamRecordImpl::MapQuality(uint8_t mapQual)
{ d_->core.qual = mapQual; return *this; }

inline PacBio::BAM::Position BamRecordImpl::MatePosition(void) const
{ return d_->core.mpos; }

inline BamRecordImpl& BamRecordImpl::MatePosition(PacBio::BAM::Position pos)
{ d_->core.mpos = pos; return *this; }

inline int32_t BamRecordImpl::MateReferenceId(void) const
{ return d_->core.mtid; }

inline BamRecordImpl& BamRecordImpl::MateReferenceId(int32_t id)
{ d_->core.mtid = id; return *this; }

inline PacBio::BAM::Position BamRecordImpl::Position(void) const
{ return d_->core.pos; }

inline BamRecordImpl& BamRecordImpl::Position(PacBio::BAM::Position pos)
{ d_->core.pos = pos; return *this; }

inline int32_t BamRecordImpl::ReferenceId(void) const
{ return d_->core.tid; }

inline BamRecordImpl& BamRecordImpl::ReferenceId(int32_t id)
{ d_->core.tid = id; return *this; }

inline bool BamRecordImpl::IsDuplicate(void) const
{ return (d_->core.flag & BamRecordImpl::DUPLICATE) != 0; }

inline BamRecordImpl& BamRecordImpl::SetDuplicate(bool ok)
{
    if (ok) d_->core.flag |=  BamRecordImpl::DUPLICATE;
    else    d_->core.flag &= ~BamRecordImpl::DUPLICATE;
    return *this;
}

inline bool BamRecordImpl::IsFailedQC(void) const
{ return (d_->core.flag & BamRecordImpl::FAILED_QC) != 0; }

inline BamRecordImpl& BamRecordImpl::SetFailedQC(bool ok)
{
    if (ok) d_->core.flag |=  BamRecordImpl::FAILED_QC;
    else    d_->core.flag &= ~BamRecordImpl::FAILED_QC;
    return *this;
}

inline bool BamRecordImpl::IsFirstMate(void) const
{ return (d_->core.flag & BamRecordImpl::MATE_1) != 0; }

inline BamRecordImpl& BamRecordImpl::SetFirstMate(bool ok)
{
    if (ok) d_->core.flag |=  BamRecordImpl::MATE_1;
    else    d_->core.flag &= ~BamRecordImpl::MATE_1;
    return *this;
}

inline bool BamRecordImpl::IsMapped(void) const
{ return (d_->core.flag & BamRecordImpl::UNMAPPED) == 0; }

inline BamRecordImpl& BamRecordImpl::SetMapped(bool ok)
{
    if (ok) d_->core.flag &= ~BamRecordImpl::UNMAPPED;
    else    d_->core.flag |=  BamRecordImpl::UNMAPPED;
    return *this;
}

inline bool BamRecordImpl::IsMateMapped(void) const
{ return (d_->core.flag & BamRecordImpl::MATE_UNMAPPED) == 0; }

inline BamRecordImpl& BamRecordImpl::SetMateMapped(bool ok)
{
    if (ok) d_->core.flag &= ~BamRecordImpl::MATE_UNMAPPED;
    else    d_->core.flag |=  BamRecordImpl::MATE_UNMAPPED;
    return *this;
}

inline bool BamRecordImpl::IsMateReverseStrand(void) const
{ return (d_->core.flag & BamRecordImpl::MATE_REVERSE_STRAND) != 0; }

inline BamRecordImpl& BamRecordImpl::SetMateReverseStrand(bool ok)
{
    if (ok) d_->core.flag |=  BamRecordImpl::MATE_REVERSE_STRAND;
    else    d_->core.flag &= ~BamRecordImpl::MATE_REVERSE_STRAND;
    return *this;
}

inline bool BamRecordImpl::IsPaired(void) const
{ return (d_->core.flag & BamRecordImpl::PAIRED) != 0; }

inline BamRecordImpl& BamRecordImpl::SetPaired(bool ok)
{
    if (ok) d_->core.flag |=  BamRecordImpl::PAIRED;
    else    d_->core.flag &= ~BamRecordImpl::PAIRED;
    return *this;
}

inline bool BamRecordImpl::IsPrimaryAlignment(void) const
{ return (d_->core.flag & BamRecordImpl::SECONDARY) == 0; }

inline BamRecordImpl& BamRecordImpl::SetPrimaryAlignment(bool ok)
{
    if (ok) d_->core.flag &= ~BamRecordImpl::SECONDARY;
    else    d_->core.flag |=  BamRecordImpl::SECONDARY;
    return *this;
}

inline bool BamRecordImpl::IsProperPair(void) const
{ return (d_->core.flag & BamRecordImpl::PROPER_PAIR) != 0; }

inline BamRecordImpl& BamRecordImpl::SetProperPair(bool ok)
{
    if (ok) d_->core.flag |=  BamRecordImpl::PROPER_PAIR;
    else    d_->core.flag &= ~BamRecordImpl::PROPER_PAIR;
    return *this;
}

inline bool BamRecordImpl::IsReverseStrand(void) const
{ return (d_->core.flag & BamRecordImpl::REVERSE_STRAND) != 0; }

inline BamRecordImpl& BamRecordImpl::SetReverseStrand(bool ok)
{
    if (ok) d_->core.flag |=  BamRecordImpl::REVERSE_STRAND;
    else    d_->core.flag &= ~BamRecordImpl::REVERSE_STRAND;
    return *this;
}

inline bool BamRecordImpl::IsSecondMate(void) const
{ return (d_->core.flag & BamRecordImpl::MATE_2) != 0; }

inline BamRecordImpl& BamRecordImpl::SetSecondMate(bool ok)
{
    if (ok) d_->core.flag |=  BamRecordImpl::MATE_2;
    else    d_->core.flag &= ~BamRecordImpl::MATE_2;
    return *this;
}

inline bool BamRecordImpl::IsSupplementaryAlignment(void) const
{ return (d_->core.flag & BamRecordImpl::SUPPLEMENTARY) != 0; }

inline BamRecordImpl& BamRecordImpl::SetSupplementaryAlignment(bool ok)
{
    if (ok) d_->core.flag |=  BamRecordImpl::SUPPLEMENTARY;
    else    d_->core.flag &= ~BamRecordImpl::SUPPLEMENTARY;
    return *this;
}

} // namespace BAM
} // namespace PacBio

#endif // BAMRECORDIMPL_H
