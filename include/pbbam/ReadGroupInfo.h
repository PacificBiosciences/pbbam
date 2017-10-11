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
/// \file ReadGroupInfo.h
/// \brief Defines the ReadGroupInfo class.
//
// Author: Derek Barnett

#ifndef READGROUPINFO_H
#define READGROUPINFO_H

#include "pbbam/Config.h"
#include "pbbam/exception/InvalidSequencingChemistryException.h"
#include <cstddef>
#include <cstdint>
#include <map>
#include <string>

namespace PacBio {
namespace BAM {

/// \brief This enum describes the base features that may be present in a read
///        group's records.
///
/// This information is stored in its description (\@RG:DS).
///
enum class BaseFeature
{
    DELETION_QV
  , DELETION_TAG
  , INSERTION_QV
  , MERGE_QV
  , SUBSTITUTION_QV
  , SUBSTITUTION_TAG
  , IPD
  , PULSE_WIDTH
  , PKMID
  , PKMEAN
  , PKMID2
  , PKMEAN2
  , LABEL
  , LABEL_QV
  , ALT_LABEL
  , ALT_LABEL_QV
  , PULSE_MERGE_QV
  , PULSE_CALL
  , PRE_PULSE_FRAMES
  , PULSE_CALL_WIDTH
  , START_FRAME
  , PULSE_EXCLUSION
};

/// \brief This enum describes the encoding types used for frame data within a
///        read group's records.
///
/// This information is stored in its description (\@RG:DS).
///
enum class FrameCodec
{
    RAW
  , V1
};

/// \brief This enum describes the experimental design of the barcodes within a
///        read group's records.
///
/// This information is stored in its description (\@RG:DS).
///
enum class BarcodeModeType
{
   NONE
 , SYMMETRIC
 , ASYMMETRIC
 , TAILED
};

/// \brief This enum describes the type of value encoded by barcode quality,
///        within a read group's records.
///
/// This information is stored in its description (\@RG:DS).
///
enum class BarcodeQualityType
{
    NONE
  , SCORE
  , PROBABILITY
};

/// \brief This enum describes the instrument type / platform model,
///        within a read group's records.
///
/// This information is stored in its description (\@RG:PM).
///
enum class PlatformModelType
{
    ASTRO
  , RS
  , SEQUEL
};

/// \brief The ReadGroupInfo class represents a read group entry (\@RG) in the
///        SAM header.
///
class PBBAM_EXPORT ReadGroupInfo
{
public:
    /// \name Conversion & Validation
    ///

    /// \brief Creates a ReadGroupInfo object from SAM-formatted text.
    ///
    /// \param[in] sam  SAM-formatted text
    /// \returns read group info object
    ///
    static ReadGroupInfo FromSam(const std::string& sam);

    /// \brief Converts a ReadGroupInfo object to its SAM-formatted text.
    ///
    /// \param[in] rg     input ReadGroupInfo object
    /// \returns SAM-formatted text (no trailing newline)
    ///
    static std::string ToSam(const ReadGroupInfo& rg);

    /// \brief Converts a read group ID (string) to its numeric value.
    ///
    /// \param[in] rgId     read group ID string
    /// \returns numeric value of ID
    ///
    static int32_t IdToInt(const std::string& rgId);

    /// \brief Converts a read group ID number to its string representation.
    ///
    /// \param[in] id     read group ID number
    /// \returns hexadecimal string representation of ID
    ///
    static std::string IntToId(const int32_t id);

    /// \returns sequencing chemistry from (bindingKig, sequencingKit,
    ///          basecallerVersion)
    ///
    static std::string SequencingChemistryFromTriple(const std::string& bindingKit,
                                                     const std::string& sequencingKit,
                                                     const std::string& basecallerVersion);

    /// \}

public:
    /// \name Constructors & Related Methods
    /// \{

    /// \brief Creates an empty read group info object.
    ReadGroupInfo();

    /// \brief Creates a read group info object with an ID.
    ///
    /// \param[in] id   string representation of read group ID
    ///
    ReadGroupInfo(std::string id);

    /// \brief Creates a read group info object from a movie name & read type.
    ///
    /// \param[in] movieName    sequencing movie name
    /// \param[in] readType     string version of record type
    ///
    /// \sa RecordType
    ///
    ReadGroupInfo(std::string movieName, std::string readType);

    /// \brief Creates a read group info object from a movie name, read type,
    ///        and platform model.
    ///
    /// \param[in] movieName    sequencing movie name
    /// \param[in] readType     string version of record type
    /// \param[in] platform     platform model type
    ///
    /// \sa RecordType
    ///
    ReadGroupInfo(std::string movieName,
                  std::string readType,
                  const PlatformModelType platform);

    ReadGroupInfo(const ReadGroupInfo&) = default;
    ReadGroupInfo(ReadGroupInfo&&) = default;
    ReadGroupInfo& operator=(const ReadGroupInfo&) = default;
    ReadGroupInfo& operator=(ReadGroupInfo&&) = default;
    ~ReadGroupInfo() = default;

    /// \}

public:
    /// \name Comparison Operators
    /// \{

    bool operator==(const ReadGroupInfo& other) const;

    /// \}

public:
    /// \name Conversion & Validation
    /// \{

    /// \returns true if read group info is valid
    ///
    /// Currently this checks to see that ReadGroupInfo::Id does not contain an
    /// empty string.
    ///
    bool IsValid() const;

    /// \brief Converts this object to its SAM-formatted text.
    ///
    /// \returns SAM-formatted text (no trailing newline)
    ///
    std::string ToSam() const;

    /// \}

public:
    /// \name Attributes
    /// \{

    /// \returns the number of barcode sequences in BarcodeFile
    ///
    /// \throws std::runtime_error if barcode data not set.
    ///         Check HasBarcodeData if this data may be absent.
    ///
    size_t BarcodeCount() const;

    /// \returns name of FASTA file containing barcode sequences
    ///
    /// \throws std::runtime_error if barcode data not set.
    ///         Check HasBarcodeData if this data may be absent.
    ///
    std::string  BarcodeFile() const;

    /// \returns MD5 hash of the contents of BarcodeFile
    ///
    /// \throws std::runtime_error if barcode data not set.
    ///         Check HasBarcodeData if this data may be absent.
    ///
    std::string BarcodeHash() const;

    /// \returns experimental design type of barcodes
    ///
    /// \throws std::runtime_error if barcode data not set.
    ///         Check HasBarcodeData if this data may be absent.
    ///
    BarcodeModeType BarcodeMode() const;

    /// \returns type of value encoded in the 'bq' tag
    ///
    /// \throws std::runtime_error if barcode data is not set.
    ///         Check HasBarcodeData if this data may be absent.
    ///
    BarcodeQualityType BarcodeQuality() const;

    /// \returns basecaller version number (e.g. "2.1")
    std::string BasecallerVersion() const;

    /// \returns tag name in use for the specified for base feature
    std::string BaseFeatureTag(const BaseFeature& feature) const;

    /// \returns binding kit part number (e.g. "100236500")
    std::string BindingKit() const;

    /// \returns true if reads are classified as spike-in controls
    bool Control() const;

    /// \returns any non-standard tags added to the \@PG entry
    ///
    /// Result map consists of {tagName => value}.
    ///
    std::map<std::string, std::string> CustomTags() const;

    /// \returns string value of \@RG:DT
    std::string Date() const;

    /// \returns string value of \@RG:FO
    std::string FlowOrder() const;

    /// \returns frame rate in Hz
    std::string FrameRateHz() const;

    /// \returns true if read group has barcode data
    bool HasBarcodeData() const;

    /// \returns true if read group has an entry for the specified base feature
    bool HasBaseFeature(const BaseFeature& feature) const;

    /// \returns string value of \@RG:ID
    std::string Id() const;

    /// \returns codec type in use for IPD
    FrameCodec IpdCodec() const;

    /// \returns string value of \@RG:KS
    std::string KeySequence() const;

    /// \returns string value of \@RG:LB
    std::string Library() const;

    /// \returns movie name (stored in \@RG:PU)
    std::string MovieName() const;

    /// \returns string value of \@RG:PL
    std::string Platform() const;

    /// \returns string value of \@RG:PM
    PlatformModelType PlatformModel() const;

    /// \returns string value of \@RG:PI
    std::string PredictedInsertSize() const;

    /// \returns string value of \@RG:PG
    std::string Programs() const;

    /// \returns codec type in use for PulseWidth
    FrameCodec PulseWidthCodec() const;

    /// \returns string value of read type
    std::string ReadType() const;

    /// \returns string value of \@RG:SM
    std::string Sample() const;

    /// \returns string value of \@RG:CN
    std::string SequencingCenter() const;

    /// \returns sequencing chemistry name
    std::string SequencingChemistry() const;

    /// \returns sequencing kit part number
    std::string SequencingKit() const;

    /// \}

public:
    /// \name Attributes
    /// \{

    /// \brief Sets read group's barcode data.
    ///
    /// Barcode fields are either absent or all must be present.
    ///
    /// \param[in] barcodeFile      barcode filename
    /// \param[in] barcodeHash      MD5 hash of barcode file
    /// \param[in] barcodeCount     number of records in barcode file
    /// \param[in] barcodeMode      experimental design of barcodes
    /// \param[in] barcodeQuality   type of barcode quality value
    ///
    /// \sa BarcodeFile \n
    ///     BarcodeHash \n
    ///     BarcodeCount \n
    ///     BarcodeMode \n
    ///     BarcodeQuality \n
    ///     ReadGroupInfo::ClearBarcodeData
    ///
    /// \returns reference to this object
    ///
    ReadGroupInfo& BarcodeData(const std::string& barcodeFile,
                               const std::string& barcodeHash,
                               size_t barcodeCount,
                               BarcodeModeType barcodeMode,
                               BarcodeQualityType barcodeQuality);

    /// \brief Sets the basecaller version number.
    ///
    /// \param[in] versionNumber   new value
    /// \returns reference to this object
    ///
    ReadGroupInfo& BasecallerVersion(const std::string& versionNumber);

    /// \brief Sets the tag to be used for a particular base feature.
    ///
    /// \param[in] feature      feature type begin updated
    /// \param[in] tag          new value
    /// \returns reference to this object
    ///
    ReadGroupInfo& BaseFeatureTag(const BaseFeature& feature,
                                  const std::string& tag);

    /// \brief Sets the binding kit part number.
    ///
    /// \param[in] kitNumber    new value
    /// \returns reference to this object
    ///
    ReadGroupInfo& BindingKit(const std::string& kitNumber);

    /// \brief Removes all barcode data from this read group.
    ///
    /// \returns reference to this read group
    ///
    ReadGroupInfo& ClearBarcodeData();

    /// \brief Removes all base features from this read group.
    ///
    /// \returns reference to this read group
    ///
    ReadGroupInfo& ClearBaseFeatures();

    /// \brief Sets whether read group's records are classifed as spike-in
    ///        controls.
    ///
    /// \param[in] ctrl     true if records are spike-in controls
    /// \returns reference to this object
    ///
    ReadGroupInfo& Control(const bool ctrl);

    /// \brief Sets a new collection of non-standard tags.
    ///
    /// Custom tag map entries should consist of {tagName => value}.
    ///
    /// \param[in] custom      new tags
    /// \returns reference to this object
    ///
    ReadGroupInfo& CustomTags(const std::map<std::string, std::string>& custom);

    /// \brief Sets the value for \@RG:DT
    ///
    /// \param[in] date      new value
    /// \returns reference to this object
    ///
    ReadGroupInfo& Date(const std::string& date);

    /// \brief Sets the value for \@RG:FO
    ///
    /// \param[in] order     new value
    /// \returns reference to this object
    ///
    ReadGroupInfo& FlowOrder(const std::string& order);

    /// \brief Sets the frame rate.
    ///
    /// \param[in] frameRateHz     string value of frame rate in Hz
    /// \returns reference to this object
    ///
    ReadGroupInfo& FrameRateHz(const std::string& frameRateHz);

    /// \brief Sets the read group's ID.
    ///
    /// \param[in] id     string value of ID
    /// \returns reference to this object
    ///
    ReadGroupInfo& Id(const std::string& id);

    /// \brief Sets the read group's ID, from movie name & read type
    ///
    /// \param[in] movieName    sequencing movie name
    /// \param[in] readType     string version of read type
    /// \returns reference to this object
    ///
    ReadGroupInfo& Id(const std::string& movieName,
                      const std::string& readType);

    /// \brief Sets the codec type used for IPD
    ///
    /// \param[in] codec    codec type
    /// \param[in] tag      IPD tag
    /// \returns reference to this object
    ///
    ReadGroupInfo& IpdCodec(const FrameCodec& codec,
                            const std::string& tag = std::string());

    /// \brief Sets the value for \@RG:KS
    ///
    /// \param[in] sequence      new value
    /// \returns reference to this object
    ///
    ReadGroupInfo& KeySequence(const std::string& sequence);

    /// \brief Sets the value for \@RG:LB
    ///
    /// \param[in] library      new value
    /// \returns reference to this object
    ///
    ReadGroupInfo& Library(const std::string& library);

    /// \brief Sets the value for movie name (stored in \@RG:PU).
    ///
    /// \param[in] movieName    new value
    /// \returns reference to this object
    ///
    ReadGroupInfo& MovieName(const std::string& movieName);

    /// \brief Sets the value for \@RG:PI
    ///
    /// \param[in] size         new value
    /// \returns reference to this object
    ///
    ReadGroupInfo& PredictedInsertSize(const std::string& size);

    /// \brief Sets the value for \@RG:PG
    ///
    /// \param[in] programs     new value
    /// \returns reference to this object
    ///
    ReadGroupInfo& Programs(const std::string& programs);

    /// \brief Sets the value for \@RG:PM
    ///
    /// \param[in] platformModel new value
    /// \returns reference to this object
    ///
    ReadGroupInfo& PlatformModel(const PlatformModelType& platform);

    /// \brief Sets the codec type used for PulseWidth
    ///
    /// \param[in] codec    codec type
    /// \param[in] tag      pulse width tag
    /// \returns reference to this object
    ///
    ReadGroupInfo& PulseWidthCodec(const FrameCodec& codec,
                                   const std::string& tag = std::string());

    /// \brief Sets the read type.
    ///
    /// \param[in] type    new value
    /// \returns reference to this object
    ///
    ReadGroupInfo& ReadType(const std::string& type);

    /// \brief Removes a particular base feature from this read group.
    ///
    /// \param[in] feature      feature to remove
    /// \returns reference to this object
    ///
    ReadGroupInfo& RemoveBaseFeature(const BaseFeature& feature);

    /// \brief Sets the value for \@RG:SM
    ///
    /// \param[in] sample       new value
    /// \returns reference to this object
    ///
    ReadGroupInfo& Sample(const std::string& sample);

    /// \brief Sets the value for \@RG:CN
    ///
    /// \param[in] center       new value
    /// \returns reference to this object
    ///
    ReadGroupInfo& SequencingCenter(const std::string& center);

    /// \brief Sets the sequencing kit part number.
    ///
    /// \param[in] kitNumber    new value
    /// \returns reference to this object
    ///
    ReadGroupInfo& SequencingKit(const std::string& kitNumber);

    /// \}

private:
    std::string id_;                    // ID * must be unique for valid SAM *
    std::string sequencingCenter_;      // CN
    std::string date_;                  // DT * (ISO-8601) *
    std::string flowOrder_;             // FO
    std::string keySequence_;           // KS
    std::string library_;               // LB
    std::string programs_;              // PG
    std::string predictedInsertSize_;   // PI
    std::string movieName_;             // PU
    std::string sample_;                // SM

    PlatformModelType platformModel_ = PlatformModelType::SEQUEL;   // PM

    // DS:<Description> components
    std::string readType_;
    std::string bindingKit_;
    std::string sequencingKit_;
    std::string basecallerVersion_;
    mutable std::string sequencingChemistry_;
    std::string frameRateHz_;
    bool        control_ = false;
    FrameCodec  ipdCodec_ = FrameCodec::V1;
    FrameCodec  pulseWidthCodec_ = FrameCodec::V1;
    bool        hasBarcodeData_ = false;
    std::string barcodeFile_;
    std::string barcodeHash_;
    size_t      barcodeCount_ = 0;
    BarcodeModeType barcodeMode_ = BarcodeModeType::NONE;
    BarcodeQualityType barcodeQuality_ = BarcodeQualityType::NONE;
    std::map<BaseFeature, std::string> features_;

    // custom attributes
    std::map<std::string, std::string> custom_; // tag => value

private:
    std::string EncodeSamDescription() const;
    void DecodeSamDescription(const std::string& description);
};

/// \brief Creates a read group ID from a movie name & read type.
///
/// \param[in] movieName    sequencing movie name
/// \param[in] readType     string version of read type
///
/// \returns hexadecimal string read group ID
///
PBBAM_EXPORT
std::string MakeReadGroupId(const std::string& movieName,
                            const std::string& readType);

} // namespace BAM
} // namespace PacBio

#include "pbbam/internal/ReadGroupInfo.inl"

#endif // READGROUPINFO_H
