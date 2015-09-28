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

#ifndef READGROUPINFO_H
#define READGROUPINFO_H

#include "pbbam/Config.h"
#include <map>
#include <string>

namespace PacBio {
namespace BAM {

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
  , LABEL
  , LABEL_QV
  , ALT_LABEL
  , ALT_LABEL_QV
  , PULSE_MERGE_QV
  , PULSE_CALL
  , PRE_PULSE_FRAMES
  , PULSE_CALL_WIDTH
};

enum class FrameCodec
{
    RAW
  , V1
};

enum class BarcodeModeType
{
   NONE
 , SYMMETRIC
 , ASYMMETRIC
};

enum class BarcodeQualityType
{
    NONE
  , SCORE
  , PROBABILITY
};

class PBBAM_EXPORT ReadGroupInfo
{
public:
    /// \name Conversion & Validation
    ///

    static ReadGroupInfo FromSam(const std::string& sam);

    static std::string ToSam(const ReadGroupInfo& rg);

    /// \}

public:
    /// \name Constructors & Related Methods
    /// \{

    ReadGroupInfo(void);
    ReadGroupInfo(const std::string& id);
    ReadGroupInfo(const std::string& movieName, const std::string& readType);
    ReadGroupInfo(const ReadGroupInfo& other);
    ReadGroupInfo(ReadGroupInfo&& other);
    ReadGroupInfo& operator=(const ReadGroupInfo& other);
    ReadGroupInfo& operator=(ReadGroupInfo&& other);
    ~ReadGroupInfo(void);

    /// \}

public:
    /// \name Attributes
    /// \{

    /// \returns the number of barcode sequences in BarcodeFile
    ///
    /// \throws std::runtime_error if barcode data not set. Check HasBarcodeData if this data may be absent.
    ///
    size_t BarcodeCount(void) const;

    /// \returns name of FASTA file containing barcode sequences
    ///
    /// \throws std::runtime_error if barcode data not set. Check HasBarcodeData if this data may be absent.
    ///
    std::string  BarcodeFile(void) const;

    /// \returns MD5 hash of the contents of BarcodeFile
    ///
    /// \throws std::runtime_error if barcode data not set. Check HasBarcodeData if this data may be absent.
    ///
    std::string BarcodeHash(void) const;

    /// \returns experimental design type of barcodes
    ///
    /// \throws std::runtime_error if barcode data not set. Check HasBarcodeData if this data may be absent.
    ///
    BarcodeModeType BarcodeMode(void) const;

    /// \returns type of value encoded in the 'bq' tag
    ///
    /// \throws std::runtime_error if barcode data is not set. Check HasBarcodeData if this data may be absent.
    ///
    BarcodeQualityType BarcodeQuality(void) const;

    /// \returns basecaller version number (e.g. "2.1")
    std::string BasecallerVersion(void) const;

    /// \returns tag name in use for the specified for base feature
    std::string BaseFeatureTag(const BaseFeature& feature) const;

    /// \returns binding kit part number (e.g. "100236500")
    std::string BindingKit(void) const;

    /// \returns true if reads are classified as spike-in controls
    bool Control(void) const;

    /// \returns additional tags not specified by either SAM/BAM or PacBio specs.
    std::map<std::string, std::string> CustomTags(void) const;

    std::string Date(void) const;

    std::string FlowOrder(void) const;

    /// \returns frame rate in Hz
    std::string FrameRateHz(void) const;

    /// \returns true if read group has barcode data
    bool HasBarcodeData(void) const;

    /// \returns true if read group has an entry for the specified base feature
    bool HasBaseFeature(const BaseFeature& feature) const;

    std::string Id(void) const;

    /// \returns codec type in use for IPD
    FrameCodec IpdCodec(void) const;

    std::string KeySequence(void) const;

    std::string Library(void) const;

    std::string MovieName(void) const;

    std::string Platform(void) const;

    std::string PredictedInsertSize(void) const;

    std::string Programs(void) const;

    /// \returns codec type in use for PulseWidth
    FrameCodec PulseWidthCodec(void) const;

    std::string ReadType(void) const;

    std::string Sample(void) const;

    std::string SequencingCenter(void) const;

    /// \returns sequencing kit part number
    std::string SequencingKit(void) const;

    /// \}

    /// \name Conversion & Validation
    /// \{

    bool IsValid(void) const;

    /// \returns SAM-formatted text representation of ReadGroupInfo
    std::string ToSam(void) const;

    /// \}

    /// \name Comparison
    /// \{

    bool operator==(const ReadGroupInfo& other) const;

    /// \}

public:
    /// \name Attributes
    /// \{

    /// Sets read group's barcode data.
    ///
    /// Barcode fields are either absent or all must be present.
    ///
    /// \param[in] barcodeFile    \sa BarcodeFile
    /// \param[in] barcodeHash    \sa BarcodeHash
    /// \param[in] barcodeCount   \sa BarcodeCount
    /// \param[in] barcodeMode    \sa BarcodeMode
    /// \param[in] barcodeQuality \sa BarcodeQuality
    ///
    /// \sa ClearBarcodeData
    ///
    /// \returns reference to this read group
    ///
    ReadGroupInfo& BarcodeData(const std::string& barcodeFile,
                               const std::string& barcodeHash,
                               size_t barcodeCount,
                               BarcodeModeType barcodeMode,
                               BarcodeQualityType barcodeQuality);

    ReadGroupInfo& BasecallerVersion(const std::string& versionNumber);

    ReadGroupInfo& BaseFeatureTag(const BaseFeature& feature,
                                  const std::string& tag);

    ReadGroupInfo& BindingKit(const std::string& kitNumber);

    /// Removes all barcode data from this read group.
    ///
    /// \returns reference to this read group
    ///
    ReadGroupInfo& ClearBarcodeData(void);

    ReadGroupInfo& Control(const bool ctrl);

    ReadGroupInfo& CustomTags(const std::map<std::string, std::string>& custom);

    ReadGroupInfo& Date(const std::string& date);

    ReadGroupInfo& FlowOrder(const std::string& order);

    ReadGroupInfo& FrameRateHz(const std::string& frameRateHz);

    ReadGroupInfo& Id(const std::string& id);

    ReadGroupInfo& Id(const std::string& movieName, const std::string& readType);

    ReadGroupInfo& IpdCodec(const FrameCodec& codec, const std::string& tag = std::string());

    ReadGroupInfo& KeySequence(const std::string& sequence);

    ReadGroupInfo& Library(const std::string& library);

    ReadGroupInfo& MovieName(const std::string& id);

    ReadGroupInfo& PredictedInsertSize(const std::string& size);

    ReadGroupInfo& Programs(const std::string& programs);

    ReadGroupInfo& PulseWidthCodec(const FrameCodec& codec, const std::string& tag = std::string());

    ReadGroupInfo& ReadType(const std::string& type);

    ReadGroupInfo& Sample(const std::string& sample);

    ReadGroupInfo& SequencingCenter(const std::string& center);

    ReadGroupInfo& SequencingKit(const std::string& kitNumber);

    /// \}

private:
    std::string id_;                     // ID * Unique ID required for valid SAM/BAM header *
    std::string sequencingCenter_;       // CN
    std::string date_;                   // DT * (ISO 8601) *
    std::string flowOrder_;              // FO
    std::string keySequence_;            // KS
    std::string library_;                // LB
    std::string programs_;               // PG
    std::string predictedInsertSize_;    // PI
    std::string movieName_;              // PU * more explicit, in place of "platform unit" *
    std::string sample_;                 // SM

    // DS:<Description> components
    std::string readType_;
    std::string bindingKit_;
    std::string sequencingKit_;
    std::string basecallerVersion_;
    std::string frameRateHz_;
    bool        control_ = false;
    FrameCodec  ipdCodec_;
    FrameCodec  pulseWidthCodec_;
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
    std::string EncodeSamDescription(void) const;
    void DecodeSamDescription(const std::string& description);
};

PBBAM_EXPORT
std::string MakeReadGroupId(const std::string& movieName,
                            const std::string& readType);

} // namespace BAM
} // namespace PacBio

#include "pbbam/internal/ReadGroupInfo.inl"

#endif // READGROUPINFO_H
