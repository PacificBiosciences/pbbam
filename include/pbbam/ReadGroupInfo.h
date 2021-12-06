#ifndef PBBAM_READGROUPINFO_H
#define PBBAM_READGROUPINFO_H

#include <pbbam/Config.h>

#include <cstddef>
#include <cstdint>

#include <map>
#include <string>
#include <utility>

#include <boost/optional.hpp>

#include <pbcopper/data/FrameCodec.h>
#include <pbcopper/data/FrameEncoders.h>
#include <pbcopper/data/Strand.h>

#include <pbbam/exception/InvalidSequencingChemistryException.h>

namespace PacBio {
namespace BAM {

/// \brief This enum describes the base features that may be present in a read
///        group's records.
///
/// This information is stored in its description (\@RG:DS).
///
enum class BaseFeature
{
    DELETION_QV,
    DELETION_TAG,
    INSERTION_QV,
    MERGE_QV,
    SUBSTITUTION_QV,
    SUBSTITUTION_TAG,
    IPD,
    PULSE_WIDTH,
    PKMID,
    PKMEAN,
    PKMID2,
    PKMEAN2,
    LABEL,
    LABEL_QV,
    ALT_LABEL,
    ALT_LABEL_QV,
    PULSE_MERGE_QV,
    PULSE_CALL,
    PRE_PULSE_FRAMES,
    PULSE_CALL_WIDTH,
    START_FRAME,
    PULSE_EXCLUSION
};

using FrameCodec PBBAM_DEPRECATED = Data::FrameCodec;

/// \brief This enum describes the experimental design of the barcodes within a
///        read group's records.
///
/// This information is stored in its description (\@RG:DS).
///
enum class BarcodeModeType
{
    NONE,
    SYMMETRIC,
    ASYMMETRIC,
    TAILED
};

/// \brief This enum describes the type of value encoded by barcode quality,
///        within a read group's records.
///
/// This information is stored in its description (\@RG:DS).
///
enum class BarcodeQualityType
{
    NONE,
    SCORE,
    PROBABILITY
};

/// \brief This enum describes the instrument type / platform model,
///        within a read group's records.
///
/// This information is stored in its description (\@RG:PM).
///
enum class PlatformModelType
{
    ASTRO,
    RS,
    SEQUEL,
    SEQUELII
};

/// \brief Aggregate to simplify ReadGroupInfo constructor.
///
struct ReadGroupInfoConfig
{
    std::string MovieName;
    std::string ReadType;
    boost::optional<PlatformModelType> Platform{};
    boost::optional<std::pair<uint16_t, uint16_t>> Barcodes{};
    boost::optional<Data::Strand> Strand{};
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

    ///
    /// \brief GetBaseId
    ///
    /// \param  id
    /// \return the hash portion only of a read group ID, with (optional)
    ///         barcode labels removed)
    ///
    /// \sa ReadGroupInfo::BaseId
    ///
    static std::string GetBaseId(const std::string& id);

    /// \brief Converts a read group ID (string) to its numeric value.
    ///
    /// \note Accepts the optional barcode-labeled IDs. These will be stripped
    ///       and number calculated from the base value.
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

    /// \brief Creates an empty read group of UNKNOWN read type.
    ReadGroupInfo();

    /// \brief Creates a read group info object with an ID.
    ///
    /// \note \p id can be a "standard" ID or contain barcode labels.
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
    ReadGroupInfo(std::string movieName, std::string readType, PlatformModelType platform);

    /// \brief Creates a read group info object with an ID.
    ///
    /// \param[in] baseId       string representation of numeric read group ID
    /// \param[in] barcodes     barcode pair for this read group
    ///
    ReadGroupInfo(std::string baseId, std::pair<uint16_t, uint16_t> barcodes);

    /// \brief Creates a read group info object from a movie name & read type.
    ///
    /// \param[in] movieName    sequencing movie name
    /// \param[in] readType     string version of record type
    /// \param[in] barcodes     barcode pair for this read group
    ///
    /// \sa RecordType
    ///
    ReadGroupInfo(std::string movieName, std::string readType,
                  std::pair<uint16_t, uint16_t> barcodes);

    /// \brief Creates a read group info object from a movie name, read type,
    ///        and platform model.
    ///
    /// \param[in] movieName    sequencing movie name
    /// \param[in] readType     string version of record type
    /// \param[in] platform     platform model type
    /// \param[in] barcodes     barcode pair for this read group
    ///
    /// \sa RecordType
    ///
    ReadGroupInfo(std::string movieName, std::string readType, PlatformModelType platform,
                  std::pair<uint16_t, uint16_t> barcodes);

    /// \brief Creates a read group info object from a ReadGroupInfoConfig
    ///
    /// \param[in] config       aggregate that contains all information to
    ///                         create a ReadGroupInfo
    ///
    /// \sa RecordType
    ///
    ReadGroupInfo(ReadGroupInfoConfig config);

    /// \}

public:
    /// \name Comparison Operators
    /// \{

    bool operator==(const ReadGroupInfo& other) const noexcept;

    /// Enable sort on RG:ID
    bool operator<(const ReadGroupInfo& other) const noexcept;

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

    /// \returns the number of barcode sequences in BarcodeFile, as stored in
    ///          the description tag (\@RG:DS)
    ///
    /// \throws std::runtime_error if barcode data not set.
    ///         Check HasBarcodeData if this data may be absent.
    ///
    size_t BarcodeCount() const;

    /// \returns name of FASTA file containing barcode sequences, as stored in
    ///          the description tag (\@RG:DS)
    ///
    /// \throws std::runtime_error if barcode data not set.
    ///         Check HasBarcodeData if this data may be absent.
    ///
    std::string BarcodeFile() const;

    /// \returns MD5 hash of the contents of BarcodeFile, as stored in the
    ///          description tag (\@RG:DS)
    ///
    /// \throws std::runtime_error if barcode data not set.
    ///         Check HasBarcodeData if this data may be absent.
    ///
    std::string BarcodeHash() const;

    /// \returns experimental design type of barcodes, as stored in the
    ///          description tag (\@RG:DS)
    ///
    /// \throws std::runtime_error if barcode data not set.
    ///         Check HasBarcodeData if this data may be absent.
    ///
    BarcodeModeType BarcodeMode() const;

    /// \returns type of value encoded by the 'bq' tag, as stored in the
    ///          description tag (\@RG:DS)
    ///
    /// \throws std::runtime_error if barcode data is not set.
    ///         Check HasBarcodeData if this data may be absent.
    ///
    BarcodeQualityType BarcodeQuality() const;

    /// \returns barcode pair stored in the read group ID (\@RG:ID)
    ///
    /// \note This does **NOT** refer to any data in the description (DS) tag.
    ///
    boost::optional<std::pair<uint16_t, uint16_t>> Barcodes() const;

    /// \returns forward barcode label stored in the read group ID (\@RG:ID)
    ///
    /// \note This does **NOT** refer to any data in the description (DS) tag.
    ///
    boost::optional<uint16_t> BarcodeForward() const;

    /// \returns reverse barcode label stored in the read group ID (\@RG:ID)
    ///
    /// \note This does **NOT** refer to any data in the description (DS) tag.
    ///
    boost::optional<uint16_t> BarcodeReverse() const;

    //// \returns string value of \@RG:BC
    std::string BarcodeSequence() const;

    /// \returns basecaller version number (e.g. "2.1")
    std::string BasecallerVersion() const;

    /// \returns tag name in use for the specified for base feature
    std::string BaseFeatureTag(BaseFeature feature) const;

    /// \returns the hash portion only of a read group ID, with (optional)
    ///          barcode labels removed)
    ///
    /// For most read groups (without barcode labels), this will be the same as
    /// ID(). However, for those read groups with barcoded-labels, this method
    /// will return the ID without those labels.
    ///
    /// Id() should be preferred over this method in most cases. This is
    /// intended for use with hash-string or integers directly.
    ///
    /// For "ID:12345678":
    ///     rg.Id()     -> "12345678"
    ///     rg.BaseId() -> "12345678"
    ///
    /// For "ID:12345678/0--0":
    ///     rg.Id()   -> "12345678/0--0";
    ///     rg.BaseId -> "12345678"
    ///
    /// \sa Id
    ///
    std::string BaseId() const;

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

    /// \returns true if the read group description (\@RG:DS) contains barcode data
    ///
    /// \note This does **NOT** refer to the optional barcode labels.
    ///
    bool HasBarcodeData() const;

    /// \returns true if read group has an entry for the specified base feature
    bool HasBaseFeature(BaseFeature feature) const;

    /// \returns full string value of \@RG:ID, whether optional barcode labels
    ///          are present
    ///
    /// This method should be perferred over BaseId() in most cases,
    /// e.g. mapping between header info.
    ///
    /// For "ID:12345678":
    ///     rg.Id()     -> "12345678"
    ///     rg.BaseId() -> "12345678"
    ///
    /// For "ID:12345678/0--0":
    ///     rg.Id()   -> "12345678/0--0";
    ///     rg.BaseId -> "12345678"
    ///
    /// \sa BaseId
    ///
    std::string Id() const;

    /// \returns codec type in use for IPD
    Data::FrameCodec IpdCodec() const;

    /// \returns codec implementation in use for IPD
    Data::FrameEncoder IpdFrameEncoder() const;

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
    Data::FrameCodec PulseWidthCodec() const;

    /// \returns codec implementation in use for PulseWidth
    Data::FrameEncoder PulseWidthFrameEncoder() const;

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

    /// \returns CCS strand
    boost::optional<Data::Strand> Strand() const;

    /// \}

public:
    /// \name Attributes
    /// \{

    /// \brief Sets barcode data for the read group's description tag.
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
    ReadGroupInfo& BarcodeData(std::string barcodeFile, std::string barcodeHash,
                               size_t barcodeCount, BarcodeModeType barcodeMode,
                               BarcodeQualityType barcodeQuality);

    /// \brief Sets the value for \@RG:BC
    ///
    /// \param[in] barcodeSequence      new value
    /// \returns reference to this object
    ///
    ReadGroupInfo& BarcodeSequence(std::string barcodeSequence);

    /// \brief Sets the basecaller version number.
    ///
    /// \param[in] versionNumber   new value
    /// \returns reference to this object
    ///
    ReadGroupInfo& BasecallerVersion(std::string versionNumber);

    /// \brief Sets the tag to be used for a particular base feature.
    ///
    /// \param[in] feature      feature type begin updated
    /// \param[in] tag          new value
    /// \returns reference to this object
    ///
    ReadGroupInfo& BaseFeatureTag(BaseFeature feature, std::string tag);

    /// \brief Sets the binding kit part number.
    ///
    /// \param[in] kitNumber    new value
    /// \returns reference to this object
    ///
    ReadGroupInfo& BindingKit(std::string kitNumber);

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
    ReadGroupInfo& Control(bool ctrl);

    /// \brief Sets a new collection of non-standard tags.
    ///
    /// Custom tag map entries should consist of {tagName => value}.
    ///
    /// \param[in] custom      new tags
    /// \returns reference to this object
    ///
    ReadGroupInfo& CustomTags(std::map<std::string, std::string> custom);

    /// \brief Sets the value for \@RG:DT
    ///
    /// \param[in] date      new value
    /// \returns reference to this object
    ///
    ReadGroupInfo& Date(std::string date);

    /// \brief Sets the value for \@RG:FO
    ///
    /// \param[in] order     new value
    /// \returns reference to this object
    ///
    ReadGroupInfo& FlowOrder(std::string order);

    /// \brief Sets the frame rate.
    ///
    /// \param[in] frameRateHz     string value of frame rate in Hz
    /// \returns reference to this object
    ///
    ReadGroupInfo& FrameRateHz(std::string frameRateHz);

    /// \brief Sets the read group's ID.
    ///
    /// \param[in] id     string value of ID
    /// \returns reference to this object
    ///
    ReadGroupInfo& Id(std::string id);

    /// \brief Sets the read group's ID, from movie name & read type
    ///
    /// \param[in] movieName    sequencing movie name
    /// \param[in] readType     string version of read type
    /// \returns reference to this object
    ///
    ReadGroupInfo& Id(const std::string& movieName, const std::string& readType);

    /// \brief Sets the codec type used for IPD
    ///
    /// \param[in] codec    codec type
    /// \param[in] tag      IPD tag
    /// \returns reference to this object
    ///
    ReadGroupInfo& IpdCodec(Data::FrameCodec codec, std::string tag = std::string());

    /// \brief Sets the codec implementation used for IPD
    ///
    /// \param[in] encoder  codec implementation
    /// \param[in] tag      IPD tag
    /// \returns reference to this object
    ///
    ReadGroupInfo& IpdFrameEncoder(Data::FrameEncoder encoder);

    /// \brief Sets the value for \@RG:KS
    ///
    /// \param[in] sequence      new value
    /// \returns reference to this object
    ///
    ReadGroupInfo& KeySequence(std::string sequence);

    /// \brief Sets the value for \@RG:LB
    ///
    /// \param[in] library      new value
    /// \returns reference to this object
    ///
    ReadGroupInfo& Library(std::string library);

    /// \brief Sets the value for movie name (stored in \@RG:PU).
    ///
    /// \param[in] movieName    new value
    /// \returns reference to this object
    ///
    ReadGroupInfo& MovieName(std::string movieName);

    /// \brief Sets the value for \@RG:PI
    ///
    /// \param[in] size         new value
    /// \returns reference to this object
    ///
    ReadGroupInfo& PredictedInsertSize(std::string size);

    /// \brief Sets the value for \@RG:PG
    ///
    /// \param[in] programs     new value
    /// \returns reference to this object
    ///
    ReadGroupInfo& Programs(std::string programs);

    /// \brief Sets the value for \@RG:PM
    ///
    /// \param[in] platformModel new value
    /// \returns reference to this object
    ///
    ReadGroupInfo& PlatformModel(PlatformModelType platform);

    /// \brief Sets the codec type used for PulseWidth
    ///
    /// \param[in] codec    codec type
    /// \param[in] tag      pulse width tag
    /// \returns reference to this object
    ///
    ReadGroupInfo& PulseWidthCodec(Data::FrameCodec codec, std::string tag = std::string());

    /// \brief Sets the codec implementation to use for PulseWidth
    ///
    /// \param[in] encoder  codec implementation
    /// \returns reference to this object
    ///
    ReadGroupInfo& PulseWidthFrameEncoder(Data::FrameEncoder encoder);

    /// \brief Sets the read type.
    ///
    /// \param[in] type    new value
    /// \returns reference to this object
    ///
    ReadGroupInfo& ReadType(std::string type);

    /// \brief Removes a particular base feature from this read group.
    ///
    /// \param[in] feature      feature to remove
    /// \returns reference to this object
    ///
    ReadGroupInfo& RemoveBaseFeature(BaseFeature feature);

    /// \brief Sets the value for \@RG:SM
    ///
    /// \param[in] sample       new value
    /// \returns reference to this object
    ///
    ReadGroupInfo& Sample(std::string sample);

    /// \brief Sets the value for \@RG:CN
    ///
    /// \param[in] center       new value
    /// \returns reference to this object
    ///
    ReadGroupInfo& SequencingCenter(std::string center);

    /// \brief Sets the sequencing kit part number.
    ///
    /// \param[in] kitNumber    new value
    /// \returns reference to this object
    ///
    ReadGroupInfo& SequencingKit(std::string kitNumber);

    /// \brief Sets the ccs strand.
    ///
    /// \param[in] strand       new value
    /// \returns reference to this object
    ///
    ReadGroupInfo& Strand(Data::Strand strand);

    /// \}

private:
    std::string id_;                   // ID * must be unique for valid SAM *
    std::string sequencingCenter_;     // CN
    std::string date_;                 // DT * (ISO-8601) *
    std::string flowOrder_;            // FO
    std::string keySequence_;          // KS
    std::string library_;              // LB
    std::string programs_;             // PG
    std::string predictedInsertSize_;  // PI
    std::string movieName_;            // PU
    std::string sample_;               // SM

    PlatformModelType platformModel_ = PlatformModelType::SEQUEL;  // PM

    // DS:<Description> components
    std::string readType_;
    std::string bindingKit_;
    std::string sequencingKit_;
    std::string basecallerVersion_;
    mutable std::string sequencingChemistry_;
    std::string frameRateHz_;
    bool control_ = false;
    Data::FrameCodec ipdCodec_ = Data::FrameCodec::V1;
    Data::FrameCodec pulseWidthCodec_ = Data::FrameCodec::V1;
    bool hasBarcodeData_ = false;
    std::string barcodeFile_;
    std::string barcodeHash_;
    size_t barcodeCount_ = 0;
    BarcodeModeType barcodeMode_ = BarcodeModeType::NONE;
    BarcodeQualityType barcodeQuality_ = BarcodeQualityType::NONE;
    std::map<BaseFeature, std::string> features_;
    boost::optional<Data::Strand> strand_;

    // (optional) barcode label handling
    boost::optional<std::pair<uint16_t, uint16_t>> barcodes_;
    std::string baseId_;

    Data::FrameEncoder ipdEncoder_ = Data::V1FrameEncoder{};
    Data::FrameEncoder pulseWidthEncoder_ = Data::V1FrameEncoder{};

    // custom attributes
    std::map<std::string, std::string> custom_;  // tag => value

private:
    std::string EncodeSamDescription() const;
    void DecodeSamDescription(const std::string& description);
    void DecodeBarcodeKey(const std::string& key, std::string value);
    void DecodeStrand(std::string value);
    std::string EncodeStrand(Data::Strand strand) const;
    void DecodeFrameCodecKey(const std::string& key, std::string value);
};

/// \brief Creates a read group ID from a movie name & read type.
///
/// \param[in] movieName    sequencing movie name
/// \param[in] readType     string version of read type
///
/// \returns hexadecimal string read group ID, e.g. "4c1bc9e4"
///
std::string MakeReadGroupId(const std::string& movieName, const std::string& readType,
                            const boost::optional<Data::Strand> strand = {});

/// \brief Creates a read group ID from a movie name, read type, and barcode string.
///
/// \param[in] movieName        sequencing movie name
/// \param[in] readType         string version of read type
/// \param[in] barcodeString    string version of barcode pair ("0--0")
///
/// \returns string containing the concatenation of the hex value with barcode
///          label "/x--y", e.g. "4c1bc9e4/0--1"
///
std::string MakeReadGroupId(const std::string& movieName, const std::string& readType,
                            const std::string& barcodeString,
                            const boost::optional<Data::Strand> strand = {});

/// \brief Creates a read group ID from a movie name, read type, and barcode IDs
///
/// \param[in] movieName    sequencing movie name
/// \param[in] readType     string version of read type
/// \param[in] barcodes     pair of barcode indices (0,0)
///
/// \returns string containing the concatenation of the hex value with barcode
///          label "/x--y", e.g. "4c1bc9e4/0--1"
///
std::string MakeReadGroupId(const std::string& movieName, const std::string& readType,
                            const std::pair<int16_t, int16_t>& barcodes,
                            const boost::optional<Data::Strand> strand = {});

/// \brief Creates a read group ID from a read group object
///
/// This convenience method detects whether barcode information is available and
/// returns the appropriate label.
///
/// \param[in] readGroup    ReadGroupInfo object
///
/// \returns string containing the concatenation of the hex value, optionally with
///          barcode label "/x--y", e.g. "4c1bc9e4" or "4c1bc9e4/0--1"
///
std::string MakeReadGroupId(const ReadGroupInfo& readGroup);

/// \brief Creates a \b LEGACY read group ID from a movie name & read type.
///
/// \warning The IDs generated by the "legacy" group of methods were incorrect, where
///          barcode information was included as part of the MD5 hash generated.
///          These are provided in case there is a need to reproduce the old behavior.
///
/// \param[in] movieName    sequencing movie name
/// \param[in] readType     string version of read type
///
/// \returns hexadecimal string read group ID
///
std::string MakeLegacyReadGroupId(const std::string& movieName, const std::string& readType);

/// \brief Creates a \b LEGACY read group ID from a movie name, read type,
///        and barcode string.
///
/// \warning The IDs generated by the "legacy" group of methods were incorrect, where
///          barcode information was included as part of the MD5 hash generated.
///          These are provided in case there is a need to reproduce the old behavior.
///
/// \param[in] movieName        sequencing movie name
/// \param[in] readType         string version of read type
/// \param[in] barcodeString    string version of barcode pair ("0--0")
///
/// \returns string containing the concatenation of the hex value with barcode label "/x--y"
///          (e.g. "4c1bc9e4/0--1")
///
std::string MakeLegacyReadGroupId(const std::string& movieName, const std::string& readType,
                                  const std::string& barcodeString);

/// \brief Creates a \b LEGACY read group ID from a movie name, read type, and
///        barcode IDs.
///
/// \warning The IDs generated by the "legacy" group of methods were incorrect, where
///          barcode information was included as part of the MD5 hash generated.
///          These are provided in case there is a need to reproduce the old behavior.
///
/// \param[in] movieName    sequencing movie name
/// \param[in] readType     string version of read type
/// \param[in] barcodes     pair of barcode indices (0,0)
///
/// \returns string containing the concatenation of the hex value with barcode label "/x--y"
///          (e.g. "4c1bc9e4/0--1")
///
std::string MakeLegacyReadGroupId(const std::string& movieName, const std::string& readType,
                                  const std::pair<int16_t, int16_t>& barcodes);

/// \brief Creates a \b LEGACY read group ID from a read group object
///
/// This convenience method detects whether barcode information is available and
/// returns the appropriate label.
///
/// \warning The IDs generated by the "legacy" group of methods were incorrect, where
///          barcode information was included as part of the MD5 hash generated.
///          These are provided in case there is a need to reproduce the old behavior.
///
/// \param[in] readGroup    ReadGroupInfo object
///
/// \returns string containing the concatenation of the hex value, optionally with
///          barcode label "/x--y", e.g. "4c1bc9e4" or "4c1bc9e4/0--1"
///
std::string MakeLegacyReadGroupId(const ReadGroupInfo& readGroup);

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_READGROUPINFO_H
